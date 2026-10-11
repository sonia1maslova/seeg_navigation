# MI 存在性森林式图

运行 `plot_pac_existence_forest.py --project-root <seeg_navigation根目录>`。

本脚本只读取 `results/pac/by_lina/formal_raw_pac` 的海马和 aMTG
`pooled_segment` 结果及已保存的 1,000 次组级 surrogate。没有新模型、
新置换或新校正；原始单侧 P 和 BH q 先独立复核再用于标注。

四个脑区沿用共享配置名称，每列六种频段组合。蓝点为真实 MI 相对
surrogate 均值的增量，灰线为中心化零分布的中央 95% 范围，**不是**
观测 MI 的置信区间。显著性不由图上范围重新判定。

两个原始纳入范围分别绘图：`all_dynamic_split` 排除终末 corridor
但保留终末 junction；`all_dynamic_split_with_terminal_C` 纳入终末
corridor。不得将这两幅图称为采用新“末端 corridor 和 junction 均排除”
规则的重分析。

输出位置：

- `results/pac/by_lina/pac_existence/pac_existence_forest_plot_data.csv`
- `figures/pac/by_lina/pac_existence/forest_pac_existence_<scope>.png`
- 同图目录的 `FDR_AND_PLOTTING_RULES.txt` 和 `forest_plot_audit.json`

原有文件不覆盖；只在明确要求替换已生成的新图时使用
`--replace-existing`。该选项不改变任何原始源文件。


## 可选零分布密度预览

`plot_pac_existence_density.py` 使用同目录 forest 脚本的只读载入与复核逻辑，显示已保存 surrogate MI 的密度、真实 MI 的位置以及原有 P/q 标记。它不运行新置换，且密度图目前是预览，不属于本次正式图库发布。需要 NumPy、pandas 和 Matplotlib。
