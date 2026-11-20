from pathlib import Path
import sys


target = Path(sys.argv[1])
text = target.read_text(encoding='utf8')
start = text.index('fn moved_round_trip(')
end = text.index('\n#[test]', start)
prefix, body, tail = text[:start], text[start:end], text[end:]
m = '    let near = |a: f64, b: f64|'
assert body.count(m) == 1
body = body.replace(m, '    let mass_first = std::env::var("MOVED_ORDER").as_deref() != Ok("recognition-first");\n' + m)
for i, name in [('f', 'faces'), ('s', 'solids')]:
    marker = f'for {i} in 0..moved.{name}.len() {{'
    assert body.count(marker) == 1
    body = body.replace(marker, f'for {i} in 0..(if mass_first {{ moved.{name}.len() }} else {{ 0 }}) {{')
marker = '    if !close(&a, &b) {'
assert body.count(marker) == 1
body = body.replace(marker, '''    let out = std::path::PathBuf::from(std::env::var_os("MOVED_OUTPUT").unwrap());
    std::fs::write(out.join("moved.json"), serde_json::to_vec(&a).unwrap()).unwrap();
    std::fs::write(out.join("read-back.json"), serde_json::to_vec(&b).unwrap()).unwrap();
''' + marker)
test = '''
#[test]
fn report_a_single_moved_slanted_bead() {
    let at = [[0.0, -1.0, 0.0, 12.5], [0.0, 0.0, -1.0, -7.25], [1.0, 0.0, 0.0, 3.0]];
    let source = fixtures().join("slanted_bead.step");
    let issues = moved_round_trip(&source, &at, true);
    let out = std::path::PathBuf::from(std::env::var_os("MOVED_OUTPUT").unwrap());
    std::fs::write(out.join("original-test-problems.json"), serde_json::to_vec(&issues).unwrap()).unwrap();
}
'''
target.write_text(prefix + body + tail + test, encoding='utf8', newline='\n')
