# 识别与证据共享分析

同时需要识别结果与面证据时，`recognise_with_evidence` 共享 STEP 读取、坐标系推断、
工作几何、候选清单和截面凹槽投影，返回两份完整文档。原有 `recognise` 与 `evidence`
接口及各自的输出保持不变。

## Python 接口

```python
import json
import quiddity

analysis = json.loads(quiddity.recognise_with_evidence("part.step"))
recognition = analysis["recognition"]
if analysis["evidence_error"] is not None:
    raise RuntimeError(analysis["evidence_error"])
evidence = analysis["evidence"]
```

参数为 `str` 文件路径或 STEP `bytes`，返回以换行结束的 JSON 字符串：

| 字段 | 内容 |
| --- | --- |
| `schema_version` | `quiddity-rust/analysis/1` |
| `recognition` | 原 `recognise(step)` 的完整返回，包括特征、指纹、坐标系及 `document` |
| `evidence` | 原证据接口的完整文档；证据被拒绝时为 `null` |
| `evidence_error` | 证据被拒绝时的原始错误文本；成功时为 `null` |

文件路径支持 `.gz`。识别的 `source_sha256` 对原始输入字节计算，包括压缩文件的压缩字节，
与原识别接口一致。读取、解压或 STEP 解析失败会抛出 `ValueError`；错误参数类型抛出
`TypeError`。输入只读取一次，两个文档来自同一份 STEP 内容。

## 坐标、类别与错误语义

识别文档的参数仍在文件坐标中。能够推断零件坐标系时，联合证据等价于
`evidence(step, framed=True)`，使用零件坐标；否则等价于 `framed=False`，使用文件坐标。
联合接口不接受独立的 `framed` 参数。需要固定文件坐标的证据消费者可继续使用原证据接口。
面、边索引和 `face_sources`、`edge_sources` 保持原始 STEP 身份，不随坐标变换重新编号。

两份文档保留各自的类别、成员、面归属、候选、拒绝记录与证明状态，不合并成同一套特征列表。
识别发布后处理可能补充面声明或撤销特征族，因此共享实现先保留发布前的候选清单，
再分别产生识别与证据。证据的组成面发现、识别指纹及两份文档的序列化仍执行。

证据生成失败不会清空识别结果，而是通过 `evidence: null` 和 `evidence_error` 明确返回。
调用者必须检查该错误，并继续遵守各文档的 `complete`、`failures` 和证明状态；
接口返回成功本身不代表所有几何结论均已获得实体证据。

正常路径只建立一次候选清单。严格模式到局部降级的重试、识别发布后处理异常所需的重新发现
均保留，不能把“一次共享调用”理解为任何输入都只进入一次发现函数。

## Rust 接口与回归

- `correspondence::recognise_with_evidence(&Part)` 返回
  `(Recognition, Result<EvidenceDocument, EvidenceViewError>)`。
- `file_api::recognise_with_evidence_file(&str)` 与
  `file_api::recognise_with_evidence_bytes(Vec<u8>)` 返回上述 JSON 封装。
- `capabilities()` 的 `functions` 列表包含新接口，`analysis` 字段声明封装版本。

`shared_analysis` Rust 测试逐字段比较原双接口与联合接口，包括变换后几何、来源、
坐标系拒绝、局部降级和故障注入；轨迹断言检查正常路径的候选发现及消费端凹槽投影各一次。
`tests/python/test_shared_analysis.py` 检查公开绑定、类型声明、路径与字节输入、摘要及错误语义。
这些测试补充原有 Rust 快照、独立几何、标签、运动不变性、数量下限及 NIST 检查，
不替代它们，不修改快照期望或放宽比较容差。
