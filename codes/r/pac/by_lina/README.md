# Lina来源的PAC分析

项目根目录：`/nfs/visu/workplace/yuhdua/projects/seeg_navigation`，对应 Windows 的 `Z:\projects\seeg_navigation`。

## 目录约定

- `codes/r/pac/by_lina/visualization/`：原结果和调整后结果的R绘图脚本。
- `codes/r/pac/by_lina/segment_effect_adjusted/`：调整终末segment排除规则的重拟合脚本。
- `results/pac/by_lina/formal_raw_pac/`：李娜原正式结果，仅迁移，不修改内容。
- `results/pac/by_lina/segment_effect_adjusted/`：本项目基于上述数据进行的R分析；包含 `exclude_terminal_corridor_and_junction` 和 `include_last` 两个平行分支，各自分为 `main_pac`、`power_control` 和 `aperiodic_control`。
- `figures/pac/by_lina/pac_segment_effect/lmm_effects/adjusted/`：当前正式结果；包含两个分析范围，每个范围下分别保存主 PAC、功率控制和 1/f 控制的热图、森林图及校正说明。
- `figures/pac/by_lina/pac_segment_effect/lmm_effects/lina_original/`：基于李娜原正式结果的对照图，保留原统计口径，不替换为调整版的六项 FDR。
- `figures/pac/by_lina/pac_segment_effect/mean_trends/lina_original/`：历史 PAC 均值趋势图，完整保留原内容。
- `figures/pac/by_lina/pac_segment_effect/mean_trends/adjusted/`：按照当前两个分析范围重新汇总的描述性均值趋势图。趋势线不是混合模型的推断结果。
- `results/pac/by_lina/pac_segment_effect/mean_trends/adjusted/`：重新汇总的被试水平和组水平均值、SEM 及核验记录，不修改原 PAC 数值。
- `codes/python/pac/by_lina/visualization/`：PAC segment 效应森林图和调整版均值趋势的计算及绘图脚本；均值计算与绘图分开。

## 2026-10-06迁移记录

|原位置|新位置|
|---|---|
|`codes/r/pac/lina_visualization`|`codes/r/pac/by_lina/visualization`|
|`results/formal_raw_pac`|`results/pac/by_lina/formal_raw_pac`|
|`figures/pac/lina_results`|`figures/pac/by_lina`|

原结果441个文件、原图132个文件、原R绘图脚本2个文件；迁移前后文件内容校验一致。检查清单位于 `results/pac/by_lina/migration_audit/`。随后仅修正本项目两份R绘图脚本的输入输出路径。原结果内部记录的历史路径不改写。

现用 Python PAC 趋势图脚本由命令参数传入CSV和输出目录，没有硬编码旧路径，不需要改动或迁移Python代码。

旧的 supplementary/legacy 绘图脚本指向 `na_li_` 的原始工作目录，与本次三处迁移无关，保持原样。

## 调整后分析

保留两个平行的分析范围。`exclude_terminal_corridor_and_junction` 按每个被试的每个试次的完整路径长度 `turns`，排除最后一个corridor的两个2.5秒子段和最后一个junction；`include_last` 则全部保留。不以筛选后剩余行的最大序号判断末段，避免把缺失数据误认作末段。

|路径条件|排除项|保留的最大segment序号|
|---|---|---|
|2次转弯|C3_first、C3_second、J2|5|
|4次转弯|C5_first、C5_second、J4|11|
|6次转弯|C7_first、C7_second、J6|17|

保留corridor子段和junction共同构成的导航顺序编号，不因信号缺失而重新编号。条件合并分析，不重新提取PAC、功率或IRASA指数。李娜原include/exclude结果始终保持不变；这里的include_last是另存的R重拟合结果。

模型：

- 主PAC：`raw_PAC ~ segment_position + (1|sub_id)`。
- 功率控制：`raw_PAC ~ segment_position + LF_power_db + HF_power_db + (1|sub_id)`。
- 非周期成分控制：`raw_PAC ~ segment_position + exponent_raw + (1|sub_id)`，分别拟合1–30 Hz和1–150 Hz。

使用Gaussian线性混合效应模型、被试随机截距和REML。为改善优化，仅在拟合内部把因变量乘1000；输出系数、标准误、区间和方差均恢复原单位。固定效应采用双侧Wald t检验，自由度为 `n - rank(X)`，与原MATLAB `fitlme`默认残差自由度口径一致；不替换为Satterthwaite检验。

R拟合对照相应分析范围的原MATLAB结果，核对系数、标准误、p值、自由度、观测数和被试数。原MATLAB的18/6/18/36项家族q值只用于实现复现核验，保存在 `audit/baseline_matlab_vs_r_full.csv`，不作为当前正式q值。

## 当前FDR校正规则：每脑区6项

按用户最新决定，BH-FDR针对 `segment_position` 固定效应的双侧原始p值，在每个脑区的6种频率组合内独立校正。6项由3个phase频段（0.5–4、4–12、12–30 Hz）与2个amplitude频段（30–70、70–150 Hz）组成。不同脑区、分析类型、指数估计范围以及include/exclude分支均不合并。

|分析|每个脑区每个集合的检验数|每个分支的集合数|模型总数|
|---|---:|---:|---:|
|主PAC：aMTG及三个海马脑区|6|4|24|
|功率控制：三个海马脑区|6|3|18|
|1/f控制：三个海马脑区分别在两个估计范围拟合|6|6|36|

每个分支共13个独立集合、78项检验。1/f的1–30 Hz和1–150 Hz分别校正：不是每脑区12项，也不是跨脑区36项。当前图仍主要显示1–150 Hz。其他固定效应保存在 `all_fixed_effects.csv`，不混入segment效应的FDR集合。

已完成的exclude模型只重新计算q值，不重复拟合，不改变原系数、标准误、p值或样本筛选；更新前的派生统计保存在 `audit/segment_effects_before_roi6_fdr.csv`。include_last从原完整数据重新进行78个R模型拟合，与相同数据的MATLAB实现核对后复用这些R拟合输出，避免重复拟合同一数据。

## 运行

```sh
Rscript codes/r/pac/by_lina/segment_effect_adjusted/run_segment_effect_adjusted.R /nfs/visu/workplace/yuhdua/projects/seeg_navigation recalculate_fdr exclude_terminal_corridor_and_junction
Rscript codes/r/pac/by_lina/segment_effect_adjusted/run_segment_effect_adjusted.R /nfs/visu/workplace/yuhdua/projects/seeg_navigation full include_last
Rscript codes/r/pac/by_lina/visualization/plot_adjusted_pac_segment_effect_heatmaps.R /nfs/visu/workplace/yuhdua/projects/seeg_navigation exclude_terminal_corridor_and_junction replace
Rscript codes/r/pac/by_lina/visualization/plot_adjusted_pac_segment_effect_heatmaps.R /nfs/visu/workplace/yuhdua/projects/seeg_navigation include_last
```

运行环境为用户指定的 `r_env`；其Rscript路径是 `/konglab/home/yuhdua/.conda/envs/r_env/bin/Rscript`。已有成功标记的模型不会重复拟合。`recalculate_fdr`明确授权仅更新已有派生结果的q值；图需要显式 `replace` 才覆盖当前分支的既定输出。

热图沿用power/encoding_lmm批准的风格：蓝白红零中心色阶、统一脑区简称、原始p值星号及BH q值短剑号；两个分支六张图共用系数色阶，频率文字靠近热图。只生成PNG、绘图数据CSV及校正说明TXT，不生成PDF/SVG。

成功标记表示计算和数值核验完成，不表示科学假设显著成立。

## 2026-10-11图目录整理

原 `figures/pac/by_lina/pac_segment_effect/heatmaps` 和 `controls` 的历史图归入 `pac_segment_effect/lmm_effects/lina_original`；原 `figures/pac/by_lina/segment_effect_adjusted` 归入 `pac_segment_effect/lmm_effects/adjusted`；原 `figures/pac/by_lina/pac_by_segment` 归入 `pac_segment_effect/mean_trends/lina_original`。原文件逐个校验后迁移，不改写历史 PNG、PDF 或 CSV。历史 manifest 内记录的当时路径保留，目录映射由根目录的整理记录说明。

现用热图脚本的默认输出位置已同步更新。混合模型结果仍保存在原 `results/pac/by_lina/segment_effect_adjusted`，李娜原结果仍保存在 `formal_raw_pac`，本次不移动或重拟合模型。

森林图直接使用各版本已有的 segment 固定效应系数和 95% 置信区间。原始版本分别核验原有 6/18/36 项 BH 集合，调整版分别核验每脑区 6 项 BH 集合。森林图星号表示原始 p 值（* < .05，** < .01，*** < .001），† 表示该版本对应的 BH q < .05。热图保留既有标注规则。

历史均值图的 18 个位置不含末 corridor，但保留末 junction。调整版纳入分支包含最长路径的 20 个位置；排除分支最多保留 17 个位置，删除各试次真实末 corridor 两个子段和末 junction，不重新编号。先在被试内平均相同 segment 位置的试次 PAC，再等权平均被试，SEM 来自被试间变异，不把触点或试次当作独立被试。

## 数值记录

每个分支的 `analysis_summary.csv` 汇总当前显著项；`fdr_manifest.csv` 与 `fdr_family_summary.csv` 明确13个六项集合；`all_segment_effects.csv` 为当前正式系数、p值和新q值。

原排除规则的78个MATLAB到R模型复现已全部通过，p和历史q值最大绝对差约为 `1.09e-7`。当前exclude的原始统计沿用已通过检查的R模型，只更改FDR分组。

`audit/old_vs_adjusted_segment_effects.csv` 保存原结果与当前结果的对照。其q值差异同时包含FDR分组改变，不能全部归因于样本排除。历史baseline大集合q仅用于实现核验；`pilot`只检查每分支一个模型，不提供完整家族FDR。
