# 李娜来源的 PAC–行为关联绘图

复用李娜正式 Stage04 中 `04C1_pac_behavior` 的24个整试次 PAC–指向误差模型。这里不重算PAC、不重新拟合模型、不修改李娜原始结果。

模型：`pointing_error ~ raw_PAC + (1|sub_id) + (1|junctions)`。PAC 是整试次、触点等权平均后的脑区 Tort MI。不是 segment 均值，也不是偏好相位向量长度。行为因变量是 pointing_error，系数为负表示更强PAC与更小指向误差相关。

## 校正

- aMTG：三相位频段 × 两振幅频段，共6项BH-FDR。
- 三个海马子区：三脑区 × 六种组合，共18项BH-FDR。
- 两个检验集合独立；只纳入04C1，不纳入04B趋势或04C2交互项。
- 整试次行为关联没有include_last/exclude_last分支，不复制成两套检验。
- 原24项校正值保留在新表的 `q_value_lina_24_tests`；`q_value` 为本次新校正值。

相位频段：0.5–4、4–12、12–30 Hz；振幅频段：30–70、70–150 Hz。

## 图

按用户提供的Power森林图样式：四脑区并排、加粗标准缩写、蓝色点和横向95%置信区间，六种相位/振幅组合由上到下排列，原始P星号和BH-FDR的†放在相应横线上方。

脑区标签从 `config/project_config.json` 读取。原始未标准化系数及区间保持不变；只调整图形显示单位，并清楚标注。星号 `* <0.05`、`** <0.01`、`*** <0.001` 表示原始P；†表示本次BH q<0.05。原95%区间不作多重比较调整。

## 路径

- 原表：`results/pac/by_lina/formal_raw_pac/04_amtg_followup_and_behavior/model_coefficients.csv`。
- 新结果：`results/pac/by_lina/pac_behavior/pac_behavior_effects_fdr_bh.csv`。
- 新图：`figures/pac/by_lina/pac_behavior/forest_pac_behavior.png`。
- R后处理：`codes/r/pac/by_lina/pac_behavior/apply_pac_behavior_fdr.R`。
- Python绘图：`codes/python/pac/by_lina/pac_behavior/plot_pac_behavior_forest.py`。
- 结果和图形目录中的TXT分别记录校正和绘图规则；JSON记录输入哈希、检验数和绘图核查。

此前的行为热图保留，不覆盖、不删除。

## 运行

```sh
env -u R_HOME /konglab/home/yuhdua/.conda/envs/r_env/bin/Rscript \
  codes/r/pac/by_lina/pac_behavior/apply_pac_behavior_fdr.R \
  /nfs/visu/workplace/yuhdua/projects/seeg_navigation

python codes/python/pac/by_lina/pac_behavior/plot_pac_behavior_forest.py \
  --project-root /nfs/visu/workplace/yuhdua/projects/seeg_navigation
```

已有校正文件时R脚本拒绝覆盖；重画本次PNG需显式提供 `--replace-existing`。
