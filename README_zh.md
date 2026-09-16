# 使用说明

本目录是重新整理的代码候选仓库，不是原始分析目录，也尚不能作为“完整一键复现包”上传。

- 代码按数据准备、RNA-CNV、YAP–Stem、Program146、功能分析、独立验证、空间验证及补充分析分组。
- `config/signatures/` 保存核对一致的 YAP17、Stem21、Program146。
- `docs/reproducibility.md` 说明当前能验证什么、缺少什么，以及如何补齐。
- `docs/spatial_reproduction.md` 说明 GSE210616 与 BSW2 两个空间数据的获取与输入准备。
- `docs/third_party_code_and_attribution.md` 区分第三方软件调用、源码修改和需要进一步追溯的实现。

代码作者与维护者：Yuting Zhang。

## 现在可以运行

在本目录运行：

```text
python tests/smoke/check_repository.py --rscript Rscript
```

这只检查代码语法、定义文件和文件校验值，不做相关性分析、不读取大型表达矩阵、不作图。通过该检查不等于全文科学结果已复现。

## 还不能跳过的步骤

1. 从工作站返回实际执行代码、参数和版本记录，确认是否发生过仅在工作站上的修改。
2. Program146 原始推导源码已恢复并归档，8 月 8 日后续版本也已匹配历史哈希。审阅 `analysis/04_program146/README.md` 中的新输入接口；完整生物学重跑及环境等价仍待验证，不能把源码恢复等同于一键复现完成。
3. 确认最终图版、补充表及尚未明确保留的扩展分析。
4. 补齐公开数据到空间分析输入对象的构建链，见 `docs/spatial_reproduction.md`。
5. 补充分析（Chen、SCAN-B、METABRIC、LUAD、CRC、HNSCC 等）的代码尚未纳入本仓库，见 `data/README.md` 的范围说明。

许可证已确定为 MIT（Copyright (c) 2026 Yuting Zhang）。论文的完整作者名单与正式引用格式将在引用定稿后补充到 `CITATION.cff`。

不要将旧 release 根目录整体上传，也不要上传内部审计目录。当前 Git 仅用于本地检查，没有创建远程仓库或推送。
