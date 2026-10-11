# 18-bin 分布聚合与复向量平均的等价性测试

代码按计算和绘图分开，默认从脚本位置识别项目根目录。

## 计算

`test_circular_vector_averaging_equivalence.py` 仅依赖 NumPy，生成模拟分布、
比较两种聚合方法，并执行缩放不变性和非法输入检查。

在项目根目录运行：

```powershell
python codes/python/pac/test/test_circular_vector_averaging_equivalence.py
```

结果保存至 `results/pac/test/circular_vector_averaging_equivalence/`：

- `simulated_distributions.csv`：6 组示例的原始值及已归一化的 18-bin 权重。
- `equivalence_results.csv`：示例与随机测试的 A/B 复向量、长度及比较结果。
- `summary.json`：随机种子、误差汇总和本次运行的脚本哈希等信息。

计算脚本拒绝覆盖非空结果目录。需要另做一轮测试时，可使用
`--project-root` 指向另一个测试根目录；不要为重画图重复计算。

## 绘图

`plot_circular_vector_averaging_equivalence.py` 依赖 NumPy 和 Matplotlib，
只读取上述两个 CSV，不生成随机数据、不重新归一化分布、不重新计算复向量。
黑色平均分布曲线由已保存的样本权重进行平均后绘制。

```powershell
python codes/python/pac/test/plot_circular_vector_averaging_equivalence.py
```

图保存至
`figures/pac/test/circular_vector_averaging_equivalence/equivalence_demo.png`。
同名 PNG 已存在时默认拒绝覆盖，确定重画时加 `--overwrite`。
还可以用 `--results-dir` 指定输入结果目录、用 `--output` 指定 PNG 路径。

## 已保存结果

2026-10-09 的拆分只调整代码组织；原先的结果表、汇总和图保持不变。
原有 `summary.json` 中的脚本哈希对应最初生成结果时的合并脚本，
不因本次代码拆分而回写历史运行信息。
