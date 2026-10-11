# 代码增量更新说明：2026-10-11

改动整理与核验于 2026-10-11 完成，代码发布于 2026-10-12。

本次整理的是上次代码发布 **`563fc1f24296311a09d95ce752411d58b232939f`** 之后的更改。更新目标为 [seeg_navigation 代码仓库](https://github.com/sonia1maslova/seeg_navigation)，不更新独立的结果图库，不更改仓库可见性，也不覆盖已有提交历史。

本次变化以可视化、目录约定和复现依赖补齐为主。**没有重新计算 PAC、PLV、Granger 或 DFI，没有重新拟合混合模型，也没有改变已有系数、原始 P 值或正式 Q 值。** PAC 新均值表是已有 PAC 指标的派生汇总，不是重新提取信号或重新计算 MI。

## 1. PLV：统一两个对齐分支的展示

- 更新 `codes/python/plv/S01_plot_plv_group_significance.py`。
- `order_aligned` 和 `terminal_aligned` 的组水平显著性图均不显示右侧海马列；这一限制只作用于绘图，不删除对应分析数据，不修改配置中的脑区列表。
- corridor 与 junction 主图不再重复嵌入图例；每种对齐模式生成一个单独、横向排列的公共图例。
- 保留五个频段的颜色、经验 P 值、原曲线和坐标口径。保留曲线已在部署时核验一致。
- 补齐此前未发布的共享 Python 配置载入器、PLV 样式与颜色 helper，以及 `config/plv_config.json`。PLV 的置换、种子、质量筛选和其他分析参数不变。

## 2. Granger：四个脑区配对的紧凑合并图

- 新增 `codes/python/granger/posthoc/group/plot_group_granger_combined_rois.py`，复用已发布的 `granger_compact_style.py`。
- 四行分别展示 aMTG 与左前海马、左后海马、左侧海马旁回及左侧后部 OFC；六列为时域 GC 和五个频段。
- 每个小图同时显示两个方向，以实心／空心方框区分，保留细的半透明连线和成对的显著触点对比例柱，不加斜线填充。
- `order_aligned` 与 `order_aligned_and_terminal_aligned` 两条分支，以及 corridor／junction，分别保存，共四张合并图。合并图是新增输出，不替换原单配对图。
- 读取现有 16 份单配对绘图 CSV 和两个参考绘图 manifest，不重新进行统计或 FDR。
- 默认输出位于 `figures/granger/formal_20261003/combined_four_rois/`；已有同名输出会触发保护，复现时建议指定新的输出目录。

## 3. DFI：指定图仅展示前三个频段

- 更新 `codes/r/dfi/bandpass/S05_plot_corridor_bandpass_dfi_group.R`。
- 对 `within_subject_condition/Angle2orientation_raw/aMTG_LAntHipp` 的 observed 图，仅显示 0.5–4、4–12 和 12–30 Hz。
- 保留原五频段版本的纵轴尺度、方向颜色、观察值、零分布均值和显著性标记；保留的 720 条绘图记录已核验一致。
- 其他脑区配对、feature 和 pointwise 图不自动改为三频段。
- 增加 `DFI_PLOT_ROI`、`DFI_PLOT_FEATURE`、`DFI_PLOT_KIND` 与 `DFI_FIGURE_ROOT`，用于限定重画对象和输出位置。限定对象的 observed 重画不会改写公共横坐标图。

## 4. PAC：目录整理、保留热图并增加森林图

- 三份 R 热图脚本更新输出目录和 provenance 位置，原模型统计与绘图内容不变。
- 所有相关图归入 `figures/pac/by_lina/pac_segment_effect/`，分为 `lmm_effects`（混合模型效应）和 `mean_trends`（描述性均值趋势）；两者均区分 `lina_original` 和 `adjusted`。
- 原始图逐文件校验后迁移；旧并行目录仅在移空后移除。原始结果 `formal_raw_pac` 和既有调整版模型结果目录不移动、不重拟合。
- 新增 `plot_pac_segment_effect_forests.py`：保留热图，同时使用已有 β、95% 置信区间、P 和 Q 绘制 12 张森林图，包含 240 个效应。置信区间不是 SEM 或贝叶斯 HPD。
- 森林图使用统一脑区简称；低、高 gamma 用绿、紫两色，公共振幅图例单独保存。星号表示原始 P 值（`* < .05`、`** < .01`、`*** < .001`），短剑号表示对应版本 BH-FDR `Q < .05`。
- **不把原始版与调整版的 Q 值混用。** 原始主分析中 aMTG 为 6 项、三个海马合计 18 项；原功率控制为 18 项；原非周期控制包含两个指数估计范围，共 36 项。调整版仍按每脑区六种频率组合分别校正，分析类型、指数范围和 include／exclude 分支各自独立。脚本对完整来源家族重新核验，但不改写正式 Q 值。

## 5. PAC：计算与绘图分离的调整版均值趋势

- 新增 `prepare_adjusted_pac_mean_trends.py`，用于派生均值表；新增 `plot_adjusted_pac_mean_trends.py`，只负责绘图。
- 输入是已有的试次、segment、脑区水平 PAC。先在同一被试、同一 segment 位置内平均试次，再在被试间等权平均；误差线为被试间 SEM。
- 新增四脑区、六种频率组合、两个分析范围共 48 张均值趋势图。corridor、junction 和总体分别显示描述性趋势；趋势线与阴影来自均值点的 OLS 拟合及 95% 均值置信带，不是新的混合模型检验。
- 原始历史均值图保留 18 个位置，仅排除末 corridor、仍包含末 junction。调整版 `include_last` 最多 20 个位置；`exclude_terminal_corridor_and_junction` 最多 17 个位置。
- 调整版按每个试次的真实路径长度排除末 corridor 的两个子段及末 junction，不是删除组曲线的最后一行；保留原导航顺序编号。
- 均值准备脚本检查调整版模型输入、筛选规则与原 PAC 一致，并对找到的历史均值表核验原始聚合口径。本次部署已核验 12 张历史均值表；如果这些表缺失，该项核验条数可能为零，不能据此宣称已完成原始均值复现。调整版运行仍需要原 PAC 与模型输入表，不是仅下载代码即可运行。

## 6. 复现依赖、边界与核验

- 新的 Python 绘图需要 NumPy、pandas、Matplotlib；PAC 均值趋势还需要 **SciPy**。R 脚本依赖按各自文件开头声明，主要包括 ggplot2、dplyr、readr、jsonlite 等。
- 设置 `SEEG_NAVIGATION_ROOT` 或使用脚本参数指定项目根目录。既有 `config/project_config.json` 保持上次发布的便携化版本和统一脑区简称，不用服务器原配置覆盖它。
- 所需 PAC 模型输入／统计表、PLV 组统计表、Granger 绘图 CSV／manifest 和 DFI 结果仍为外部输入；本次不上传真实数据表、信号、个体名单、拟合对象、置换数组、图像或凭证。
- 上次发布时已经修正的 MI-density 说明和非周期 β 水平绘图说明保留远端版本，不用服务器上尚未同步的旧 README 回退它们。
- 仓库中已有的 `codes/Figure*.py`、历史复现数据和其他旧文件保留；它们在当前服务器代码树中没有同名文件，并不意味着应当删除。
- 本次来源文件有 SHA256 和 Git blob 校验记录；Python 语法及共享配置载入、R 语法解析在发布前检查。先前部署完成的图形／数值核验作为验证依据，不为发布重复运行正式分析。

参见 [PAC 目录与复现指南](PAC_SEGMENT_EFFECT_DIRECTORY_GUIDE.md) 和 [本次来源快照](source_snapshot_20261011_update.json)。上次的 [发布说明](RECENT_ANALYSIS_AND_FIGURES_20261011.md) 保留为历史记录。
