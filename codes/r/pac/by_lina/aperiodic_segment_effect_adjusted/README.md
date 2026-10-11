# 李娜来源的非周期性谱指数 segment 趋势补充分析

本目录研究 **1/f 指数本身**随 segment 序号的变化，不是把指数作为协变量的 PAC 控制分析。

## 当前分析

- 脑区：aMTG、左前海马、左后海马、右前海马。
- 主要指标：已保存的 1–150 Hz IRASA 非周期性谱指数；不重算 IRASA。
- 模型：`exponent_raw ~ segment_position + condition + (1|sub_id)`。
- `condition` 为分类固定效应，水平来自项目配置：`turns_num2`、`turns_num4`、`turns_num6`；以 2 次转弯为参照。用原试次编号交叉核对 condition 与 turns。
- corridor 的两个 2.5 s 子段和 junction 均纳入，按原导航顺序连续编号，筛选后不重新编号。不额外加入 condition 交互项或 segment_type。
- 采用 REML 和 `lme4::lmer`，被试随机截距。双侧系数 t 检验自由度沿用此前 R 对照 MATLAB 的 `n − rank(X)` 口径，保存相应 95% t 置信区间。
- 拟合阶段保留原始 P 值，不自动做 FDR。按用户随后确认的方案，使用独立后处理脚本：aMTG 每条分支只有一个目标检验，保留原始 P；左前、左后、右前三个海马子区在每条分支内共同做 BH-FDR（每次 3 个检验）。两分支分别校正，不把 condition 或截距系数纳入。

两条平行分支：

1. `include_last`：保留全部末 corridor 和末 junction。
2. `exclude_terminal_corridor_and_junction`：每个被试的每个试次，根据完整路径的 turns 排除 `C(turns+1)_first`、`C(turns+1)_second` 和 `J(turns)`。不能用信号筛选后剩余最大编号判断末段。

## 输入与聚合

输入是李娜已有的独立指数表 `pipeline3_exponent_trial_level_dataset.csv` 的原样快照，位于新结果根目录的 `inputs/lina_exponent_trial_level_dataset.csv`。原路径和 SHA256 保存在 `analysis_manifest.json`。

每一行对应被试 × 脑区 × 试次 × 2.5 s 子段 × 指数拟合范围。源表中的触点指数已等权平均为脑区水平，模型不再次做触点聚合。不使用 PAC 六个频段组合中的重复指数行。读取后仅选择四脑区的 1–150 Hz 范围，并记录非有限指标的剔除情况。

`exponent = −slope`：segment 系数为正代表指数增大，对应 log–log 谱斜率变得更负。

## 输出

结果根目录：`results/pac/by_lina/aperiodic_segment_effect_adjusted/`。

- `all_segment_effects.csv`：两分支共八个 segment 系数及其原始 P、标准误和区间。
- `all_segment_effects_fdr_bh.csv`：保留原始系数和 P，新增三海马 BH 校正值、检验数和显著性依据。aMTG 的 `p_fdr_bh` 留空，`p_for_significance` 使用原始 P；海马使用 BH 校正 P。
- `FDR_CORRECTION_NOTES.txt`、`fdr_analysis_manifest.json` 和 `_FDR_SUCCESS`：记录后处理政策、结果及原始文件未改动的 SHA256。原 `analysis_manifest.json` 中的“不做 FDR”仅描述拟合阶段，保留不改写。
- `all_fixed_effects.csv`：全部固定效应，包括两个 condition 对比。
- `model_diagnostics.csv`：观测数、被试数、条件数、模型收敛、奇异性和方差。
- 每分支：分析数据、四个 RDS 模型、相应统计和按 condition 的数据覆盖表。
- `audit/`：行筛选记录、各条件删除数、原始输入和既有结果的前后 SHA256。
- `_SUCCESS`：八个模型完成、基础数值校验通过、受保护文件保持不变；不代表研究假设显著。

原李娜结果和之前 `segment_effect_adjusted` 的 PAC 模型不修改。此分析不启动后台监控。后续按用户确认生成紧凑系数图，两分支各一张 PNG，位于 `figures/pac/by_lina/aperiodic_segment_effect_adjusted/`。

## 系数图

绘图代码：`codes/python/pac/by_lina/aperiodic_segment_effect_adjusted/plot_aperiodic_segment_effect_forest.py`，只读取 R 模型结果，不拟合模型。

采用同一个坐标轴，纵轴为四脑区统一缩写，横轴为原模型 β；蓝色点表示估计值，水平误差线表示95%置信区间。横轴使用科学计数法标注 `×10⁻³`，估计值本身不变。标记紧贴各点及其区间：星号均按原始 P（`* <0.05`、`** <0.01`、`*** <0.001`）；`†` 单独表示三个海马在相应分支内共同 BH-FDR 后 q<0.05。aMTG不校正、不加†。不打印底部显著性说明，规则单独保存为 TXT。区间为原模型未作多重比较调整的95%区间。两分支使用相同横轴范围以便比较。

两张 PNG 分别放入 `include_last/` 和 `exclude_terminal_corridor_and_junction/` 子目录。无 PDF/SVG；绘图清单、所用星号、输入 SHA256 和只读核查位于图形目录的 `plot_audit.json`。

## 运行

```sh
env -u R_HOME /konglab/home/yuhdua/.conda/envs/r_env/bin/Rscript \
  codes/r/pac/by_lina/aperiodic_segment_effect_adjusted/run_aperiodic_segment_effect_adjusted.R \
  /nfs/visu/workplace/yuhdua/projects/seeg_navigation
```

使用用户指定的 `r_env`。已有成功结果或非空模型分支时脚本会停止，不自动覆盖、删除或重跑。

完成拟合后，单独运行校正（不会重新拟合）：

```sh
env -u R_HOME /konglab/home/yuhdua/.conda/envs/r_env/bin/Rscript \
  codes/r/pac/by_lina/aperiodic_segment_effect_adjusted/apply_aperiodic_segment_fdr.R \
  /nfs/visu/workplace/yuhdua/projects/seeg_navigation
```

已有校正结果时后处理脚本同样拒绝覆盖。

用具有 numpy、pandas、matplotlib 的 Python 运行绘图：

```sh
python codes/python/pac/by_lina/aperiodic_segment_effect_adjusted/plot_aperiodic_segment_effect_forest.py \
  --project-root /nfs/visu/workplace/yuhdua/projects/seeg_navigation
```

默认不覆盖已存在的两张图；显式 `--replace-existing` 仅用于重绘这两张已命名的 PNG。
