# PAC segment 效应：目录与复现指南

调整版作为当前正式分析，李娜原始版本作为对照；两版使用同一来源的 PAC 指标，但保留各自筛选和 FDR 口径。

## 图和结果的组织

所有路径均相对于项目根目录：

```text
figures/pac/by_lina/pac_segment_effect/
├── lmm_effects/
│   ├── lina_original/
│   │   ├── include_last/
│   │   ├── exclude_last/
│   │   └── provenance/
│   └── adjusted/
│       ├── include_last/
│       └── exclude_terminal_corridor_and_junction/
└── mean_trends/
    ├── lina_original/
    └── adjusted/
        ├── include_last/
        └── exclude_terminal_corridor_and_junction/
```

`lmm_effects` 每个分析范围下分别保存 `main_pac`、`power_control` 和 `aperiodic_control_1_150Hz`，每类下面保留 `heatmaps` 并增加 `forest_plots`。模型结果本身仍在 `results/pac/by_lina/segment_effect_adjusted/`，李娜原始结果仍在 `results/pac/by_lina/formal_raw_pac/`。

新派生均值表位于 `results/pac/by_lina/pac_segment_effect/mean_trends/adjusted/<branch>/`，含被试水平、组水平、SEM、覆盖统计与核验记录；不是新的 MI 计算或模型拟合。

历史 PNG、PDF、CSV 和 manifest 迁移时保持原内容。历史 manifest 中的旧路径是生成时记录，不追溯改写。已部署图目录的 `DIRECTORY_MIGRATION_MANIFEST.json` 记录新旧位置映射。

## 脚本与运行顺序

新的 Python 脚本均位于 `codes/python/pac/by_lina/visualization/`，共享脑区简称来自 `config/project_config.json`。先设置项目根目录：

```sh
export SEEG_NAVIGATION_ROOT=/path/to/seeg_navigation
```

在原始与调整版模型统计表均已准备好的条件下，生成森林图：

```sh
python codes/python/pac/by_lina/visualization/plot_pac_segment_effect_forests.py \
  --project-root "$SEEG_NAVIGATION_ROOT" \
  --output-root "$SEEG_NAVIGATION_ROOT/figures/pac/by_lina/reproduced_segment_forests"
```

默认森林图根目录就是正式 `pac_segment_effect`。建议复现时指定新的目录，避免覆盖；只有显式 `--replace-existing` 才允许覆盖既有森林图输出。该脚本不会重新拟合模型或修改统计输入。

均值计算和绘图分两步：

```sh
python codes/python/pac/by_lina/visualization/prepare_adjusted_pac_mean_trends.py \
  --project-root "$SEEG_NAVIGATION_ROOT"

python codes/python/pac/by_lina/visualization/plot_adjusted_pac_mean_trends.py \
  --project-root "$SEEG_NAVIGATION_ROOT" \
  --output-root "$SEEG_NAVIGATION_ROOT/figures/pac/by_lina/reproduced_mean_trends"
```

运行前需提供原完整 ROI PAC 表和两个分支的调整版模型输入。另需提供已迁移的历史 PAC 均值表，才能完成原始均值的对照核验：脚本对找到的表逐一检查，但历史表全缺失时这项核验会记录零条，不会自动报错。本次部署已核验 12 张历史表；复现时应查看 manifest 的实际核验条数，不把零条核验理解为对照已通过。

现用 R 热图脚本仍位于 `codes/r/pac/by_lina/visualization/`：原始主效应、原始两种控制，以及调整版效应分别生成热图；本次只更新其输出目录，不改变模型或校正公式。相关模型运行与统计规则详见 `codes/r/pac/by_lina/README.md`。

## 读图与统计口径

- 森林图为既有模型 β 和 95% 置信区间，低／高 gamma 用绿／紫，公共图例单独保存。星号为原始 P，短剑号为该版本的 BH Q。
- 原始 Q 使用来源表中的完整 6／18／36 项校正家族；调整版每个脑区内六项，不跨分析类型、指数范围或 include／exclude 分支合并。
- 描述性均值先被试内平均试次，再被试等权平均；误差线为被试间 SEM。趋势线的 95% 阴影为组均值点 OLS 拟合的均值置信带，不表示混合模型置信区间或显著性。
- 原历史均值图有 18 个位置，仍保留末 junction；调整 include 分支最多 20 个位置，exclude 分支最多 17 个位置。按真实试次路径判断末段并保持原编号，不根据当前数据表最后一行判断。

需要 Python 3.10 或更新版本，以及 NumPy、pandas、Matplotlib、SciPy。这里没有新增环境锁文件，也没有在发布时重新运行正式模型。新图默认保存 PNG 和必要数据／核验表，不额外保存 PDF 或 SVG；迁移的历史 PDF 原样保留。
