"""复用资源仪器，平衡比较原基线、当前 PR 和细窄面修复候选。"""

import hashlib
import json
from pathlib import Path
import sys

import arm_resources_run as shared


def orders():
    stages = ("base", "pr1", "candidate")
    forward = [stages[i:] + stages[:i] for i in range(3)]
    return (forward + [tuple(reversed(order)) for order in forward]) * 2


if __name__ == "__main__":
    plan = json.loads(Path(sys.argv[sys.argv.index("--plan") + 1]).read_text())
    assert plan["adapter_sha256"] == hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    shared.STAGES = ("base", "pr1", "candidate")
    shared.orders = orders
    shared.main()
