# Rankine 源一阶水动力数值审计报告

生成时间：2026-08-30 05:26:17（Asia/Shanghai）

## 审计边界与原则

本轮只审计并修正单球一阶 Rankine 源 BEM、网格/缓存调用链和基础验证。二阶漂移未参与参数调节，也未用于本轮通过判据。

所有 CRESTU/WAMIT 数值均为原始物理输出；没有使用经验 scale factor、强制拟合、频率相关修正或绘图层替换。

## Forensic 网格调用链

1. `Run_SingleSphere_Convergence` 解析每层 compute/reload config，并通过 `run_frequency_domain_case` 进入同一个生产求解入口。
2. `build_bmf_domain` 读取真实 body BMF，提取 body waterline，并分别生成/读取 free-surface、bottom、outer-boundary panels；外域参考水线与 body 水线已解耦。
3. 每个频率重新构造活动域；panel centroid 同时作为 collocation point，全部 panel vertices/normal 作为 Rankine source-panel 几何。source/collocation 均做 SHA-256 并写入审计表。
4. `assemble_rankine_matrix` 对每个 collocation–source panel 对调用 `rankine_panel_integrals`，装配 D/S、自项 0.5、自由面 Robin、海底 Neumann 与外边界辐射条件。
5. radiation RHS 为 `S*(i*omega*nj)`，diffraction RHS 为 `S*dphi_I/dn`；二者均由 `solve_complex_system` 解同一个实际影响矩阵，并保存相对线性残差。
6. 逐频率 audit 的 source/collocation count、unknown count、matrix rows/columns 必须一致；本轮全部通过。

## 已确认的 bug

1. 旧 Coarse/Medium/Fine 的物面确实不同（108/300/588 panels），但自由面同时为 120/400/1064、海底为 228/700/1652、外边界为 72/120/168 panels。因此旧试验混合了 body discretization 与 Rankine outer-domain error。
2. 旧结果复用仅核对 schema、频率、浪向和漂移字段，没有核对实际 BMF 内容；网格改变后仍可能错误复用结果。
3. 旧势缓存使用面板数与若干数组求和作为 geometry signature，未包含逐频率自由面/海底/外域的抗碰撞哈希，也未完整包含 DOF、problem type、solver options 和代码版本。
4. 旧一阶流程把 pressure-integrated radiation damping 计算出来后，又用 energy/Haskind damping 覆盖最终 `B`；这会把验证量误当主求解量，并可直接污染 RAO。
5. 旧 `Run_SingleSphere_Convergence` 的源码路径和 WAMIT 根目录解析少了一层 `Source Code`/版本目录，入口依赖于外部 MATLAB path 状态。

## 已修复的 bug

1. Body-mesh 组现在强制使用 3/5/7 cube-face divisions，并验证实际 panels、unique vertices、characteristic size 与 SHA-256。
2. Body-mesh 组使用 Fine 球体水线作为共同 outer-domain reference，逐频率自由面、海底、外域哈希在三层之间完全一致。
3. Outer-domain 组固定同一个 Fine body mesh，只细化自由面/海底/外边界；body hash 在三层之间完全一致。
4. 势缓存升级为 schema 5；基础键与逐频率键均包含四类 mesh hash、频率、浪向、DOF、problem type、solver options、对称组和数值源码版本。
5. 每个频率强制输出并保存 source/collocation counts、BEM unknowns、matrix size、component hashes、cache key 与 hit/miss。
6. `A_raw/B_raw`、`A_reciprocal/B_reciprocal` 与独立 `damping_energy` 同时保存；最终 RAO 使用 pressure damping 的互易投影，不再被 energy damping 覆盖。
7. Fine clean recompute 与 cache reload 逐字段执行机器精度对账。

## 原验证中无效的结论

- 旧三层曲线不能用于声称 body mesh convergence，因为三类 outer-domain boundaries 同时改变。
- 旧结果文件不能证明当时配置对应的网格真实参与求解，因为结果级 cache compatibility 未包含 mesh fingerprint。
- 旧 `B33` 与 RAO 的趋势不能单独归因于物面网格，因为最终 damping 来自另一条 energy/Haskind 路径且外域同时变化。
- 互易投影后的零残差不能证明原始离散系统满足 reciprocity；必须查看本报告保存的 raw residual。

## 新的实际求解网格

| Study | Level | Body panels | Vertices | Waterline | FS | Bottom | Outer | BEM unknowns |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BodyMesh | Coarse | 108 | 121 | 24 | 280 | 868 | 112 | 1368 |
| BodyMesh | Medium | 300 | 321 | 40 | 280 | 868 | 112 | 1560 |
| BodyMesh | Fine | 588 | 617 | 56 | 280 | 868 | 112 | 1848 |
| OuterDomain | Coarse | 588 | 617 | 56 | 112 | 700 | 56 | 1456 |
| OuterDomain | Medium | 588 | 617 | 56 | 280 | 868 | 112 | 1848 |
| OuterDomain | Fine | 588 | 617 | 56 | 560 | 1148 | 224 | 2520 |

逐频率完整记录见 `Case_Wave/Case1_SingleSphere_Convergence/SingleSphere_Case_Audit.csv`。

## 新的真实网格收敛结果

### 固定 outer-domain 的 body mesh convergence

| Quantity | RMS Medium error | RMS Fine error | median abs(MF)/abs(CM) | Fine closer fraction | Fine overall closer |
|---|---:|---:|---:|---:|---:|
| A33 | 1.152668e-02 | 1.148002e-02 | 0.280362 | 0.938 | 1 |
| B33 | 1.893905e-02 | 1.477139e-02 | 0.204747 | 0.562 | 1 |
| F3 | 2.149017e-02 | 1.237858e-02 | 0.276752 | 1.000 | 1 |
| RAO3 | 4.393367e-01 | 2.502595e-01 | 0.276171 | 1.000 | 1 |

### 固定 Fine body 的 outer-domain convergence

| Quantity | RMS Medium error | RMS Fine error | median abs(MF)/abs(CM) | Fine closer fraction | Fine overall closer |
|---|---:|---:|---:|---:|---:|
| A33 | 2.265243e+01 | 1.660231e+01 | 0.432486 | 0.938 | 1 |
| B33 | 2.005166e+03 | 4.306451e+01 | 14.418629 | 0.938 | 1 |
| F3 | 4.045152e+01 | 2.466903e+01 | 0.623865 | 0.938 | 1 |
| RAO3 | 5.082818e+01 | 9.781254e+01 | 0.607155 | 0.875 | 0 |

频率曲线、相邻网格差值与 quantity/error–h 图：

- `Case_Wave/Case1_SingleSphere_Convergence/BodyMesh_Frequency_Curves.png`
- `Case_Wave/Case1_SingleSphere_Convergence/BodyMesh_Grid_Differences.png`
- `Case_Wave/Case1_SingleSphere_Convergence/BodyMesh_h_Convergence.png`
- `Case_Wave/Case1_SingleSphere_Convergence/OuterDomain_Frequency_Curves.png`
- `Case_Wave/Case1_SingleSphere_Convergence/OuterDomain_Grid_Differences.png`
- `Case_Wave/Case1_SingleSphere_Convergence/OuterDomain_h_Convergence.png`

这里的零网格尺度值是三点二次外推，只用于内部离散趋势检查，不作为 WAMIT 拟合目标。

## 基础 verification tests

| Level | Area rel. error | Volume rel. error | Buoyancy-center error (m) | Surface identity residual | Outward normals | Pass |
|---|---:|---:|---:|---:|---:|---:|
| Coarse | 1.484983e-02 | 2.953549e-02 | 1.677328e-02 | 2.953549e-02 | 1.000 | 1 |
| Medium | 5.384295e-03 | 1.074690e-02 | 6.085942e-03 | 1.074690e-02 | 1.000 | 1 |
| Fine | 2.752500e-03 | 5.499325e-03 | 3.111826e-03 | 5.499325e-03 | 1.000 | 1 |

- Raw A reciprocity residual 最大值：1.140996e-03。
- Raw B reciprocity residual 最大值：2.353983e-03。
- Fine outer 高段 pressure/energy B33 residual 中位数：42.610%。
- 辐射/绕射线性方程最大相对残差：6.047517e-14 / 1.012037e-13。

### Fine clean recompute / cache reload

| Field | Max absolute difference | Machine tolerance | Pass |
|---|---:|---:|---:|
| A_raw | 0.000000e+00 | 1.862645e-09 | 1 |
| B_raw | 0.000000e+00 | 3.725290e-09 | 1 |
| A_reciprocal | 0.000000e+00 | 1.862645e-09 | 1 |
| B_reciprocal | 0.000000e+00 | 3.725290e-09 | 1 |
| damping_energy | 0.000000e+00 | 3.725290e-09 | 1 |
| excitation | 0.000000e+00 | 7.450581e-09 | 1 |
| rao_complex | 0.000000e+00 | 9.094947e-13 | 1 |

## 仍未解释的 WAMIT 偏差

| Quantity | Full-grid median error | High-frequency median error |
|---|---:|---:|
| A33 | 19.581% | 26.783% |
| B33 | 77.422% | 72.921% |
| F3 | 8.880% | 53.397% |
| RAO3 | 18.223% | 58.615% |

上述偏差来自 Fine body + Fine outer 的直接无缩放比较，详见 `Case_Wave/Case1_SingleSphere_Convergence/Fine_FirstOrder_WAMIT_Comparison.csv`。

## 高频 B33/F3/RAO 退化定位

在 omega >= 1.5 rad/s，Medium→Fine 的 body sensitivity（A33/B33/F3/RAO3）中位数分别为 0.213% / 0.563% / 0.683% / 1.197%。
同一区间的 outer-domain sensitivity 分别为 6.504% / 96.051% / 54.459% / 60.050%。
pressure 与 energy B33 的高段中位差为 42.610%，因此二者不一致仍是独立、未通过的物理验证项；它不再被用来覆盖 pressure damping。

## 下一步最可能的问题位置

1. 若 pressure/energy residual 显著大于离散 sensitivity，优先审计 Haskind/energy prefactor、有限水深 depth factor、相位约定和远场积分归一化。
2. 若 outer sensitivity 大于 body sensitivity，优先检查逐频率 sponge 半径、自由面 Robin 条件和 far-boundary 辐射条件；不要继续细化物面来掩盖。
3. 若 F3 在两类网格都稳定但仍偏离 WAMIT，优先检查 diffraction Neumann 边界符号、incident/diffraction symmetry decomposition 与波幅归一化。
4. 若 A/B raw reciprocity residual 不随 h 下降，优先检查 source/collocation 自项、面元法向以及 Rankine panel integral 的 D/S 符号，而不是再做互易投影。
5. 单体一阶稳定后，再进入双体交叉辐射耦合；near/far-field mean drift 继续保持后置。

## 可复现产物

- `Case_Wave/Case1_SingleSphere_Convergence/SingleSphere_Convergence_Summary.mat`
- `Case_Wave/Case1_SingleSphere_Convergence/FirstOrder_Verification.csv`
- `Case_Wave/Case1_SingleSphere_Convergence/Fine_Cache_Clean_Reload_Verification.csv`
- `Case_Wave/Case1_SingleSphere_Convergence/Sphere_Mesh_Verification.csv`
- `Case_Wave/Case1_SingleSphere_Convergence/BodyMesh_Convergence_Long.csv`
- `Case_Wave/Case1_SingleSphere_Convergence/OuterDomain_Convergence_Long.csv`
