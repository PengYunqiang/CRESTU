# CRESTU-1F Phase 3.2 预声明网格与判据（v1，冻结）

冻结时间：2026-09-01（Asia/Shanghai）

适用 checkpoint：`7f2d391d1b42366fdb04e1b493439dab438584eb`

状态：在读取任何 Phase 3.2 新网格数值结果之前声明。

本文件 v1 一旦产生任一新网格数值结果即不可覆盖。若 v1 后来证明不适用，只能保留 v1、新建带版本号的文件、说明理由、重跑全部受影响算例，并在 reviewer gate 中记录。

## 1. 固定物理与数值公式

以下量在 L1/L2/L3 间固定：直径 10 m 的解析半球、水深 50 m、`rho=1025 kg/m^3`、`g=9.80665 m/s^2`、质量/重心/惯性/静水恢复、DOF 顺序 `[surge,sway,heave,roll,pitch,yaw]`、参考点、全域无对称缩减、source formulation、radiation/diffraction/traction convention、outer ABC、incident-wave convention、WAMIT reference 数值和每个频率的物理 outer extent。

每个频率先由冻结 Phase 2.3 基线 `q_top=2.0`、`q_bottom=1.5` 得到 `R_top,base` 和 `R_bottom,base`，再固定

```text
R_top,Phase3.2    = 1.5 * R_top,base
R_bottom,Phase3.2 = 1.5 * R_bottom,base
q_level(omega)    = (R_Phase3.2(omega) - r_inner) / lambda(omega)
```

同一频率的三个等级使用完全相同的 `R_top,Phase3.2` 和 `R_bottom,Phase3.2`。网格级别不得改变物理外边界。

## 2. 协同非嵌套网格 family

```ini
NESTED_MESH = NO
MESH_FAMILY = coordinated_non_nested_L1_L2_L3_v1
```

| 控制量 | MESH_L1 | MESH_L2 | MESH_L3 | 作用区域 |
|---|---:|---:|---:|---|
| body cube-face division | 3 | 5 | 7 | body surface |
| FS radial counts `[near,sponge]` | `[3,2]` | `[4,3]` | `[6,4]` | free surface / annulus |
| bottom radial counts `[near,sponge]` | `[3,2]` | `[4,3]` | `[6,4]` | bottom / annulus |
| configured wall vertical count | 2 | 3 | 4 | FARFIELD slant target |
| common FARFIELD theta count | 56 | 80 | 108 | wall angular discretization |
| FS outer theta count | 56 | 80 | 108 | FS outer transition |
| bottom outer theta count | 44 | 64 | 84 | bottom outer transition |
| mesh transition radius | 20 m | 20 m | 20 m | fixed geometry control |
| sponge start radius | 20 m | 20 m | 20 m | fixed physical/numerical control |

所有 theta count 均为 4 的倍数。三个等级必须由 production generator 重新生成，不复制原仓库三个 modified `OuterDomain` BMF。

控制空间的联合标称尺寸定义为七个倒数控制量的几何平均：

```text
eta = geomean([1/n_body, 1/n_FS_radial_total, 1/n_bottom_radial_total,
               1/n_wall_vertical, 1/n_wall_theta,
               1/n_FS_theta, 1/n_bottom_theta])
eta_L1 = 0.0901298343593916
eta_L2 = 0.0614843623387160
eta_L3 = 0.0448408251673150
r21 = eta_L1 / eta_L2 = 1.4658984972938
r32 = eta_L2 / eta_L3 = 1.3711692884620
```

`eta` 仅用于协同 family 的 Richardson 横坐标；它不是 total-panel-count 替代品。只有当 body/FS/bottom/wall/transition 的实际 characteristic spacing 全部随等级下降、各 component hash 和 operator hash 均变化时，才允许使用它。

## 3. 实际网格尺度、geometry 和 isolation

每个 mesh/frequency 至少保存 body、FS、bottom、wall 的 `mean(sqrt(area))`、`median(sqrt(area))`，transition-region 的相同量，以及仅作成本摘要的 global equivalent `h=sqrt(total area/total panels)`。不能只用 total panel count 声称收敛。

每级必须证明：四类 panel count、实际 characteristic spacing、mesh/operator hash 均变化；解析半球半径、水深、outer extent、质量/静水和 production physics 未变化。任一未隔离项使 `MESH_CONTROL_ISOLATION=FAIL`。

Geometry gate 使用以下预声明阈值：

- 所有顶点、中心、法向、面积有限，面积严格为正；有效最短边 `>1e-10 m`；duplicate panel count 为 0。
- 解析半球顶点半径误差 `<=1e-8 m`，body 法向外指；FS 法向 `+z`；bottom 法向 `+z`；FARFIELD 法向具有负径向分量并满足冻结 frustum convention。
- body/FS waterline最大间隙 `<=1e-8 m`；FS/wall 和 bottom/wall 环闭合误差 `<=1e-8 m`。
- 忽略三角形 BMF 重复的第 4 顶点后，closed-domain 非零 edge multiplicity 必须为 2；Euler/manifold、duplicate vertex/panel 和 analytic self-intersection audit 必须显式报告。
- 质量阈值沿用现有 `check_mesh_quality` 语义：body aspect ratio `<=4.0`，其余近场 `<=4.5`，sponge/transition/wall `<=12.0`；skewness `<=55 deg`，corner angle `[35,145] deg`，warp `<=8 deg`。单一 component poor-panel fraction `>5%` 判 FAIL。

## 4. 求解、algebraic 和 provenance gate

每个 mesh/frequency 使用唯一临时目录和 case name，配置 `IPOTEN=1`，不得读取旧 cache。结果提取完成后只保留正式 CSV/figure/report；临时 BMF/MAT/cache 不提交。

```ini
minimum_raw_rcond = 1e-12
minimum_scaled_rcond = 1e-12
maximum_relative_linear_or_algebraic_residual = 1e-10
maximum_abs_outer_ntheta = 1e-10
maximum_outer_column_reconstruction_residual = 1e-14
```

任一 production function、cfg、BMF、cache、result 或 report 解析到原仓库，立即令 `RUNTIME_PROVENANCE_GATE=FAIL` 并停止昂贵 sweep。任何 NaN/Inf 或 algebraic gate failure 令对应 case FAIL，并停止扩大 sweep 以保留最小复现。

## 5. Observable、scale floor 与 complex 比较

保存 `A33`、`B33`、complex `F3`、complex `RAO3`、complex dynamic denominator `D3`、raw/scaled rcond 和四类 residual。采用半球半径 `R=5 m`：

```text
A33_floor = 1e-6 * rho * R^3                         [kg]
B33_floor = 1e-6 * rho * omega * R^3                 [kg/s]
F3_floor  = 1e-6 * rho * g * R^2                     [N/m]
RAO3_floor = 1e-6                                    [m/m]
D3_floor  = 1e-6 * max(|-omega^2(M33+A33)|, |C33|)   [N/m]
```

相邻 complex quantity 使用

```text
relative_complex_difference = |Q_fine-Q_coarse| / max(|Q_fine|, Q_floor)
```

并同时保存 complex difference、absolute difference、magnitude difference、wrapped/unwrapped phase difference和所用 floor。只有 `|Q_fine|>100*Q_floor` 且 unwrap 连续时，相位才参与 gate；低幅值相位仅报告。

## 6. 相邻敏感性与 mesh status

非共振 observable 的 L2→L3 gate：

- `A33/B33/F3 relative difference <=2%`：PASS；`>2%` 且 `<=10%`：PARTIAL；`>10%`：FAIL。
- `RAO3 relative difference <=5%`：PASS；`>5%` 且 `<=15%`：PARTIAL；`>15%`：FAIL。
- `omega=1.4 rad/s` 的单点 RAO 不用于单独判定共振收敛，改由第 8 节 resonance tracking 判定。

`SENTINEL_CONVERGENCE` 和 `FORMAL_CONVERGENCE` 取其覆盖范围内的最坏状态。原始 non-monotonic 和 FAIL 行必须保留。

## 7. Observed order、Richardson 与 GCI

只对 A33、B33、Re/Im(F3)、Re/Im(RAO3) 等明确实标量尝试 generalized unequal-ratio Richardson。要求：三个有限值、相邻差同号且非零、无符号跨零、`eta_L1>eta_L2>eta_L3`、在 `p in [0.1,10]` 内有唯一根，并且该 component 的实际 h 单调下降。

计算：

```text
(Q1-Q2)/(Q2-Q3) = (eta1^p-eta2^p)/(eta2^p-eta3^p)
Q_ext = Q3 + (Q3-Q2)/(r32^p-1)
GCI32 = 1.25 * |(Q3-Q2)/Q3| / (r32^p-1)
```

fine-to-extrapolated error和 GCI 均保留。只有 `0.5<=p<=6`、fine adjacent sensitivity `<=5%` 且 generalized asymptotic ratio位于 `[0.8,1.2]` 时，该 observable 才记 `ASYMPTOTIC_RANGE=DEMONSTRATED`。否则写 `NOT_DEMONSTRATED`，并分类为 oscillatory convergence、pre-asymptotic、resonance shift、mesh-family inconsistency 或 unresolved。禁止直接对 complex number 套用实数 GCI。

## 8. 共振跟踪

每个 mesh level 固定运行 `omega=1.30:0.02:1.50 rad/s`，覆盖冻结 1.4 rad/s 峰值两侧。使用同样的 clean/isolated production solve。比较 peak frequency、peak `|RAO3|`、complex RAO trajectory、`D3` minimum frequency、minimum `|D3|`、峰值附近 A33/B33、damping contribution、amplitude continuity和 unwrapped phase continuity。

L2→L3 resonance gate：peak frequency shift `<=0.02 rad/s`、peak magnitude difference `<=10%`、minimum `|D3|` difference `<=10%`、峰值 B33 difference `<=10%`、相邻扫描点 amplitude step `<=50%` 且有效幅值区 unwrapped phase step `<=120 deg`。超出 2 倍上述阈值为 FAIL；介于一倍和两倍之间为 PARTIAL。

## 9. 执行顺序与资源降级

先运行 sentinel `{0.5,0.6,1.0,1.4,1.5,2.0} rad/s × 3 levels`。无 implementation-level blocker 后运行 resonance scan，再运行 `0.5:0.1:2.0 rad/s × 3 levels`。

只有 sentinel 实测外推显示完整 L3 formal sweep 将超过 2 小时、任一 case 超过 20 分钟，或可靠峰值内存估计超过 24 GiB，才允许降级为三等级 sentinel加 L2/L3 formal。必须保存证据；两级 formal 不计算/声称完整 Richardson/GCI，且不放宽 gate。

## 10. Production mesh recommendation 与 Task 7 前置条件

推荐 L3 必须同时满足 runtime provenance、mesh control isolation、geometry、algebraic、sentinel、formal 和 resonance gate 均非 FAIL，且所有非共振 L2→L3 observable 满足 PASS 阈值；若 asymptotic range 不是普遍 demonstrated，只能标记 `PASS_WITH_LIMITATIONS`，不能写 production ready。

进入定量 WAMIT 归因必须同时满足：

```makefile
MESH_GEOMETRY_GATE = PASS
ALGEBRAIC_GATE = PASS
SENTINEL_CONVERGENCE != FAIL
FORMAL_CONVERGENCE != FAIL
READY_FOR_QUANTITATIVE_WAMIT_AUDIT = YES
WAMIT_REFERENCE_PROVENANCE = PASS
```

外部 WAMIT 原始文件公开再分发状态若未确认，只允许提交 hash、derived CSV、metadata、figure 和 report，不提交原始文件。

## 11. 明确禁止的事后调整

禁止 smoothing、scale factor、empirical correction、threshold relaxation、删/藏 FAIL 频点、根据 WAMIT 调物理参数、把 non-monotonic 结果改写成单调、改变 radiation/diffraction RHS、traction、A/B extraction、ABC、source orientation、RAO 公式、WAMIT 数值、统一换行、renormalize、修改 `.gitattributes`、复用不匹配 cache、silent fallback、复制原仓库 modified OuterDomain BMF，以及把 PARTIAL 写成 production ready。
