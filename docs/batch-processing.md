# 全量处理与 pathway 文件夹

批量目录为 `data/reactome-98-bulk/`，与原有 demo、原下载目录及旧链检查点分开。数据管线只使用确定性 Python、冻结配置和英文模板。

```text
reactome-98-bulk/
  manifest.json
  batch_status.json
  processing_config.json
  pathway_hierarchy.json
  pathway_index.jsonl
  index.sqlite
  catalog/
    raw/                         # 原始层级、物种、pathway 清单和补查响应
    receipts.json                # 每个响应的 URL、SHA-256、版本
    pathway_membership.jsonl     # 包含子 pathway 的 reaction 归属
    status.json
  global/
    snapshot.json
    raw/
    downloads.sqlite
    download_status.json
    normalized/
    graph/
    paths/
      paths.sqlite
      complete_chains.jsonl
      cycle_stopped_chains.jsonl
      status.json
    chains.txt                   # 所有已导出全局链的文本，空行分隔
    audit/
  by_species/
    Bos_taurus__taxon-9913/
      Glutathione_synthesis_and_recycling__R-BTA-174403/
        manifest.json
        reactions.jsonl
        edges.jsonl
        chain_refs.jsonl         # 全局 chain ID 与当前 pathway 的成员反应 ID
```

## 下载和续跑

所有命令从仓库根目录执行。当前已有批量目录，无须再次执行 seed。

```bash
python3 -m reactome_pipeline.bulk status --dataset data/reactome-98-bulk
python3 -u -m reactome_pipeline.bulk fetch --dataset data/reactome-98-bulk
```

fetch 校验已有原始 JSON 和回执，只下载缺失记录。原始版本固定为 98，查询参数使用冻结的 enhanced/v2 nested 配置。发生网络错误时保存进度；恢复网络后重新执行相同命令。不要切换到不同 release 的在线服务继续同一快照。

吞吐由 `--workers` 和 `--interval` 控制：`--interval` 是所有线程共用的请求启动间隔，硬上限为 1/interval 次/秒（默认 0.2 秒即 5 次/秒）。实测 12 个线程 + `--interval 0.12` 约 8 条/秒；提高并发前先确认服务端状态。520–524（Cloudflare 源站错误）与其他 5xx 按指数退避重试，连续失败 8 次后保存进度并停止。

```bash
python3 -u -m reactome_pipeline.bulk fetch --dataset data/reactome-98-bulk --workers 12 --interval 0.12
```

pathway 层级是独立下载任务：

```bash
python3 -u -m reactome_pipeline.catalog --dataset data/reactome-98-bulk
```

如果该任务仍在运行，目录锁会阻止第二个任务同时写入。当前启动记录在 `runs/bulk-jobs.json`，日志分别为 `runs/bulk-fetch.log` 和 `runs/bulk-catalog.log`。日志中的历史错误保留；判断进度时同时查看进程和状态文件。

目录索引先冻结全部 Pathway 清单，获取各主浏览物种的完整 event hierarchy，再逐个查询未被这些树覆盖的 Pathway/TopLevelPathway。成员关系由源字段计算并向父 pathway 传播，共享子 pathway 保留全部父节点；无法解析的成员会停止索引构建。没有 pathway 归属的 reaction 保留在全局数据，并按记录的物种归入 `unassigned` 浏览组。

新建批量目录时，可复用已有下载：

```bash
python3 -m reactome_pipeline.bulk seed --source data/reactome-98 --dataset data/reactome-98-bulk
```

seed 保留快照、原始文件和下载回执，不复制旧的 normalized、graph 或 paths。原始文件优先建立硬链接，跨文件系统时复制；后续写入采用原子替换。

## 全量标准化、构图和枚举

reaction 下载和 catalog 都完成后执行：

```bash
python3 -u -m reactome_pipeline.bulk process --dataset data/reactome-98-bulk \
  --pathway R-BTA-174403
```

依次执行标准化、原始字段审计、全局构图、图审计、pathway 索引和链枚举。原始字段审计不通过时停止。`--pathway` 只选择自动导出的浏览文件夹，不限制下载、构图或链枚举。

枚举优先处理小型弱连通分量，然后处理大型分量。一个工作批次最多运行 1,000,000 次 DFS 操作、300 秒遍历时间或新增 10,000 条链；达到任一工作预算后保存检查点、审计结果并更新索引，接着继续下一批。这些预算不是链长度限制。时间预算在 DFS 操作之间检查，不包括随后审计和导出的时间。

只运行一批检查结果：

```bash
python3 -u -m reactome_pipeline.bulk process --dataset data/reactome-98-bulk \
  --pathway R-BTA-174403 --rounds 1
```

不加 `--rounds` 会持续推进，直到枚举完成、用户中断、出错或磁盘剩余空间低于默认 10 GiB 预留值。Ctrl-C 在枚举阶段会保存进度。`processing_config.json` 冻结代码校验值、模板版本、处理顺序和工作预算；后续复用同一批次时应保持这些参数。`--rounds` 可变化。

`batch_status.json` 显示阶段；`global/paths/progress.json` 显示正在遍历的进度；`global/paths/status.json` 和根目录 `manifest.json` 记录已经导出的批次状态。`enumeration_complete: false` 表示链仍未穷举完成。所有链结束后的完整性状态与 reaction 下载完整性分开记录。

每次枚举批次退出时同步生成 `global/chains.txt`，先列 sink 链，再列遇环停止的链，各组保持枚举输出顺序。文本直接读取结构化记录的 `text`，以流式方式写入并原子替换，不把全部链载入内存。尚未开始枚举时不生成该文件；枚举未完成时，它包含截至当前批次已导出的链。

## 浏览和按需导出

`pathway_index.jsonl` 列出全部 pathway 的名称、目标目录、成员数与当前可用链数。全局链只保存一份；SQLite 索引增量记录 chain ID 与 pathway 的关系。

```bash
python3 -m reactome_pipeline.organize export --dataset data/reactome-98-bulk \
  --pathway R-BTA-174403
```

命令会更新索引，然后生成四个文件：manifest、reactions、edges 和 chain_refs。它与处理任务使用同一目录锁，应在处理批次退出后运行；持续处理时可以提前通过 process 的 `--pathway` 参数指定自动更新的文件夹。

文件范围：

- `reactions.jsonl`：当前 pathway 及其子 pathway 的成员反应。
- `edges.jsonl`：内部边和跨边界边，用 `view_scope`、`source_in_pathway`、`target_in_pathway` 标记。
- `chain_refs.jsonl`：每条记录包含 `chain_id` 和 `member_reaction_ids`；不保存完整链结构或文本。
- `manifest.json`：来源、物种、文件校验值、成员覆盖情况、索引进度及全局枚举是否完成。

按引用读取时保留完整链，外部反应不截断，可从全局标准化记录按 ID 查找。一个全局 chain ID 可出现在多个 pathway 文件夹中；汇总链数量时按 ID 去重。物种目录用于浏览上下文，不参与全局实体匹配。图保持完整；未导出的 pathway 仍可通过索引选择导出。

## 查看全文与显式导出副本

直接从全局数据库读取某个 pathway 的完整链，输出到终端：

```bash
python3 -m reactome_pipeline.organize show --dataset data/reactome-98-bulk \
  --pathway R-BTA-174403
```

查看命令不生成链全文文件，索引状态输出到 stderr，链文本输出到 stdout。只有明确需要分享独立文件时才运行：

```bash
python3 -m reactome_pipeline.organize export-full --dataset data/reactome-98-bulk \
  --pathway R-BTA-174403 --out runs/glutathione-chain-copy
```

`export-full` 的目标必须在数据集目录之外，输出 `chains.jsonl` 和 `chains.txt`。默认的 pathway 导出不调用它。升级旧 pathway 文件夹时，仅移除与旧 manifest 校验值一致的自动生成链副本；修改过或没有回执的文件会保留并报错。

## 离线检验目录结构

```bash
python3 scripts/create_browsing_demo.py --dataset runs/browsing-demo
python3 -m unittest discover -v
```

该示例仅使用冻结的完整谷胱甘肽 pathway，包含 7 个反应、11 条边和 10 条链，不访问网络。全库处理不以这个 pathway 为范围。
