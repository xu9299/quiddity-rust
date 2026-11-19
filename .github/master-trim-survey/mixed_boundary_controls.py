"""Check the identity of both survey drivers against analytic STEP controls."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys
import tempfile

baseline, candidate, expected_file = map(Path, sys.argv[1:])
text = Path("crates/haecceity/tests/fixtures/master-trim-square.step").read_text()
plane, count = re.subn(r"(?m)^#108=.*?;$", "#108=PLANE('',#11);", text)
assert count == 1
results = {}
with tempfile.TemporaryDirectory() as workspace:
    path = Path(workspace) / "one-face.step"
    for height, chosen_area in [("0.04", 1 - 0.04 / 3), ("2.4", 1.)]:
        bowed = plane.replace("(0.5,0.0004)", f"(0.5,{height})")
        mixed = bowed.replace("#142=SURFACE_CURVE('',#135,(#141),.PCURVE_S1.)", "#142=SURFACE_CURVE('',#135,(#141),.CURVE_3D.)")
        foreign = bowed.replace("#117=PCURVE('',#108", "#197=PLANE('',#11);\n#117=PCURVE('',#197")
        assert len({bowed, mixed, foreign}) == 3
        for label, control in [("local", bowed), ("mixed", mixed), ("other_carrier", foreign)]:
            raw = control.encode()
            path.write_bytes(raw)
            digest = hashlib.sha256(raw).hexdigest()
            values = [json.loads(subprocess.check_output([exe, path, digest]))["faces"][0]["area"] for exe in (baseline, candidate)]
            assert abs(values[0] - 1.) < 1e-12 and abs(values[1] - chosen_area) < 1e-12, (label, height, values)
            results[f"{label}_{height}"] = {"input_sha256": digest, "reference_area_mm2": chosen_area, "baseline_and_candidate": values}
expected_file.write_text(json.dumps(results, indent=2) + "\n")
