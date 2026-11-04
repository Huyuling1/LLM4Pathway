# LLM4Pathway

将 Reactome 全库反应整理为英文反应文本、带类型的反应图和完整候选链，为后续 SFT 数据构建提供输入。运行过程使用 Python 标准库、固定模板和图算法，不调用 LLM。

## 运行

需要 Python 3.11+，使用 Linux 或 macOS。在仓库根目录运行：

```bash
python3 -m reactome_pipeline run --out data/reactome
```

该命令依次执行全库发现、下载、规范化、构图、图审计、链枚举和结果审计。首次发现覆盖所有 `ReactionLikeEvent` 子类、全部物种、推断反应及疾病相关反应。下载任务按清单中的稳定顺序调度，默认 4 个线程，全局请求启动间隔为 0.2 秒。

本次 release 98 的全库清单为 96,065 个反应：[清单摘要](reports/reactome-98-inventory-summary.json)，[完整清单及快照配置](reports/reactome-98-inventory.json.gz)。摘要记录各反应子类的数量和清单指纹。

再次运行同一命令会复用已校验的下载文件并继续链枚举。查看进度：

```bash
python3 -m reactome_pipeline status --out data/reactome
```

也可以分阶段执行：

```bash
python3 -m reactome_pipeline discover  --out data/reactome
python3 -m reactome_pipeline fetch     --out data/reactome
python3 -m reactome_pipeline normalize --out data/reactome
python3 -m reactome_pipeline graph     --out data/reactome
python3 -m reactome_pipeline audit     --out data/reactome
python3 -m reactome_pipeline paths     --out data/reactome
python3 -m reactome_pipeline audit     --out data/reactome
```

`Ctrl+C` 会停止下载调度，或在当前 DFS 操作结束后保存枚举位置。重启时使用同一个输出目录。写入命令对该目录加锁；`status` 可在任务运行时读取进度。

## 本地示例

六个固定的合成反应展示物质连接、催化连接、正负调控、分支、环及货币分子排除：

```bash
python3 examples/create_example.py --out runs/example
python3 -m reactome_pipeline normalize --out runs/example
python3 -m reactome_pipeline graph     --out runs/example
python3 -m reactome_pipeline paths     --out runs/example
python3 -m reactome_pipeline audit     --out runs/example
```

该示例生成 7 条带类型的证据边、6 条到达图中末端的链、2 条遇环停止的记录。预生成结果在 [examples/expected/synthetic](examples/expected/synthetic)。

仓库还包含 release 98 的两个真实反应快照：`PIP3 activates mTORC2` 和 `mTORC2 phosphorylates SGK1`。离线处理命令：

```bash
python3 examples/create_example.py --real --out runs/real-example
python3 -m reactome_pipeline normalize --out runs/real-example
python3 -m reactome_pipeline graph     --out runs/real-example
python3 -m reactome_pipeline paths     --out runs/real-example
python3 -m reactome_pipeline audit     --out runs/real-example
```

前一个反应的产物 `R-HSA-9980228` 与后一个反应的催化实体相同，形成一条 `catalyst_dependency` 边。ATP、ADP 保留在文本中。结果在 [examples/expected/real](examples/expected/real)。这些小例用于检查输出；全库范围由 `discover` 获取。

## 代码架构

```text
reactome_pipeline/
├── cli.py          命令行、阶段调度、完成状态
├── fetch.py        API、全库清单、版本检查、缓存和重试
├── normalize.py    JSOG 引用解析、计量关系和反应字段提取
├── currency.py     冻结名单及其指纹
├── graph.py        SQLite 精确连接、证据导出、强连通分量
├── paths.py        最大简单路径枚举、遇环记录、事务检查点
├── render.py       英文反应与链文本模板
├── audit.py        覆盖、拓扑、排除统计、precedingEvent 对照
├── locking.py      每个快照目录的写锁
├── common.py       规范 JSON、SHA-256、原子文件写入
└── config/
    └── currency.v1.json
```

## 数据来源和规范化

全库清单来自 `data/schema/ReactionLikeEvent/min`，数量与 `data/schema/ReactionLikeEvent/count` 核对。API 按名称分页，同名记录可能跨页重复；发现阶段改变分页边界，对 `dbId` 取并集，直至数量与服务端一致。最终清单按 `stId, dbId` 排序并冻结。

每个反应单独请求：

```text
GET /ContentService/data/query/enhanced/v2/{dbId}
    ?view=nested
    &includeRef=true
    &fetchIncomingRelationships=true
    &summariseReferenceEntity=true
```

原始响应压缩保存为 `raw/reactions/{dbId}.json.gz`，收据记录解压后内容的 SHA-256 和请求 URL。下载前后检查 Reactome release；重启时检查已有文件的内容指纹。API release 变化时，应新建快照目录。

规范化先建立 `@id → 对象` 索引，再解析 `@ref`、`inputs/outputs` 的 `element` 包装、化学计量、催化及调控关系。只将根反应的参与者提取到角色列表中。实体的名称、区室、修饰状态使用源字段；无 `stId` 的参与者写入问题清单，不参与稳定 ID 匹配。

每个规范化反应包含：

```text
id / db_id / st_id / st_id_version / schema_class / display_name
inputs[] / outputs[] / catalysts[] / regulators[]
species[] / species_name / is_inferred / inferred_from[]
pathways[] / preceding_events[] / negative_preceding_events[]
source {file, sha256, url, query, reactome_release}
text
```

参与者保存 `entity_id`、`display_name`、`compartments`、`chebi_ids`、`stoichiometry` 和 `source_pointer`。调控者另外保存原始调控类与正负号。

## 反应文本

固定格式为：

```text
inputs → outputs; catalyst: ...; positive regulators: ...; negative regulators: ...
```

例如真实快照中的反应：

```text
ATP [cytosol] + SGK1 [cytosol] → ADP [cytosol] + p-S422-SGK1 [cytosol]; catalyst: mTORC2:PIP3 [plasma membrane]
```

参与者按稳定标识及源指针排序。源数据缺少输入、输出或催化实体时写 `[not recorded]`；缺失化学计量在 JSON 中保存为 `null`。

## 连边规则

对每一个上游输出实体，在全库下游角色中匹配相同的 physical-entity `stId`：

| 上游角色 | 下游角色 | 边类型 |
|---|---|---|
| output | input | `material_flow` |
| output | catalystActivity.physicalEntity | `catalyst_dependency` |
| output | 正调控 regulator | `positive_regulation` |
| output | 负调控 regulator | `negative_regulation` |
| output | 类型未明确的 regulator | `unknown_regulation` |

`Requirement` 归入正调控，并保留原始类。复合物和集合以其自身 `stId` 匹配。不同 compartment、修饰状态或复合物身份维持为不同实体。Pathway 和物种归属作为元数据保留，不作为连边过滤条件；跨 pathway、跨物种的候选数量在 audit 中报告。

一条证据边对应 `(source, target, kind, entity_id)`，包含内容指纹及两端角色的 JSON 指针。同一反应对可有多种证据；枚举时只走一条邻接关系，链中附上全部证据。自环保留在图中。

### 货币分子

[currency.v1.json](reactome_pipeline/config/currency.v1.json) 冻结水、质子、ATP/ADP/AMP、无机磷酸及焦磷酸、NAD/NADH、NADP/NADPH 的明确 ChEBI 标识，包括已列出的带电形式。也可填写特定 Reactome 实体 ID。

按参与者 `referenceEntity` 的标识匹配，不按名称或复合物成员匹配。排除覆盖所有下游角色，反应文本保持完整。构图时复制配置并保存指纹，`excluded_entities.jsonl` 记录被排除的实体和对应候选边数。修改名单后使用 `--currency path/to/config.json` 重新构图；已有路径检查点必须与图指纹一致。

## 完整链与环

枚举单位是**两端都无法再加入未出现 reaction 的有向简单路径**。这里的“最大”表示不能继续延长。孤立 reaction 保存为单节点链。

算法从无上游节点及所有环中节点开始，按稳定 ID 顺序深度优先遍历。当前链用集合检查是否重复 reaction；走到无法继续延伸的末端后，再检查起点能否接入未出现的前序 reaction。可继续向前延伸的片段不导出。

```text
R1 → R2 → R4           保存为 sink
      ↓ ↑
      R3

R1 → R2 → R3          保存为 cycle_blocked，记录 R3 → R2
```

不会生成包含第二次 R2 的链，例如 `R1 → R2 → R3 → R2 → R4`。纯环按不同起点保存有限序列。每条链保存 reaction ID 序列、连接证据和模板文本；调控连接行明确显示正负调控。

默认穷举，没有路径长度上限或随机采样。枚举状态与链记录在同一 SQLite 事务中提交，每 1,000 次 DFS 操作或约 15 秒保存检查点。可分批运行：

```bash
python3 -m reactome_pipeline paths --out data/reactome --stop-after 10000
python3 -m reactome_pipeline paths --out data/reactome --stop-after 10000
```

每次增加最多 10,000 条新记录，再次执行继续同一遍历。`--stop-after` 控制暂停位置，未结束的枚举保持 `enumeration_complete: false`。

## 输出

```text
snapshot.json                     冻结清单、release、API 及查询配置
raw/inventory/                    清单的原始分页响应
raw/reactions/*.json.gz           每个反应的原始 nested JSON
downloads.sqlite                 下载收据
download_status.json             下载覆盖率
download_errors.jsonl             本次运行失败记录
normalized/reactions.jsonl        规范化反应及英文文本
normalized/issues.jsonl           缺失标识、名称及未解析引用
normalized/status.json           规范化覆盖率
graph/graph.sqlite               全量图和角色关系
graph/edges.jsonl                 带类型、实体和源指针的证据边
graph/currency.frozen.json        本次构图使用的名单
graph/excluded_entities.jsonl     货币分子排除统计
graph/status.json                图规模、指纹和输入覆盖状态
paths/paths.sqlite               链及事务检查点
paths/progress.json              当前 DFS 位置
paths/complete_chains.jsonl       到达无后继节点的最大简单链
paths/cycle_stopped_chains.jsonl  因重复 reaction 而停止的最大简单链
paths/status.json                枚举状态和数量
audit/report.json                覆盖率、拓扑和字段对照
audit/cyclic_components.jsonl    环所在的强连通分量
audit/preceding_only.jsonl        仅 precedingEvent 支持的连接
audit/matched_rejected_connections.jsonl  匹配到人工否决的连接
```

`precedingEvent` 与 `negativePrecedingEvent` 用于对照审计，不改变上述匹配规则。图中无后继表示当前图的边界；`sink` 与 `cycle_blocked` 分开统计。

全库结果完成需要同时满足：下载清单覆盖完整、规范化覆盖完整、`paths/status.json` 中 `enumeration_complete` 与 `complete` 都为 `true`，且审计中的 `paths_match_current_graph` 为 `true`。`--allow-partial` 可处理当前已下载部分，覆盖状态会传递到图和链。退出码 `0` 表示该阶段完成，`2` 表示覆盖或枚举仍未完成，`1` 表示运行错误。

原始全库数据和运行中的 SQLite 不提交 Git；仓库保存代码、配置、小型输入快照、预生成示例和全库清单摘要。

## 验证

```bash
python3 -m unittest discover -v
```

测试覆盖 JSOG 前向引用、三类连接及调控符号、货币分子、计量文本、缺失字段、下载恢复、内容校验、release 变化、环、孤立节点、路径恢复和部分覆盖状态。路径算法另与独立的排列穷举程序核对全部 64 种三节点无自环有向图，并加入自环情况。

## 来源

- [Reactome Content Service](https://reactome.org/dev/content-service)
- [Reactome 数据模型](https://reactome.org/documentation/data-model)
- [参与分子的字段和复合物结构](https://reactome.org/dev/graph-database/extract-participating-molecules)
- [Reactome Data Model Glossary](https://download.reactome.org/documentation/DataModelGlossary_V90.pdf)
- [ReactionT5](https://link.springer.com/article/10.1186/s13321-025-01075-4)：角色化反应文本的参考；本项目使用 Reactome 名称与固定模板。
