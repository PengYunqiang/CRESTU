# CRESTU-1F Phase 3.2 预声明网格与判据（v2，冻结）

冻结时间：2026-09-01（在查看任何 v2 网格结果之前）

适用 checkpoint：`7f2d391d1b42366fdb04e1b493439dab438584eb`

## 1. 版本关系和触发证据

v1 文件 `Phase3_2_Predeclared_Gates.md` 和全部 v1 原始 FAIL CSV 永久保留。v1 的 18 个 mesh-only case 证明 analytic geometry、outer extent、四区域 panel count、component hash、edge multiplicity、Euler/manifold、orientation 和 interface closure 均正确，但 mesh family 未通过 v1 已声明的 quality/isolation gate：

- explicit outer theta count 与各级 legacy ring count 不匹配，count-transition 三角形出现最小角约 4.78 deg、最大 aspect ratio 122.20；
- L1 body 的 maximum warp 为 8.99 deg，poor-panel fraction 为 5.56%；
- 0.5/0.6 rad/s 的 FS/bottom median characteristic spacing 非单调。

证据见 `Phase3_2_V1_Mesh_Manifest.csv`、`Phase3_2_V1_Mesh_Quality.csv`、`Phase3_2_V1_Mesh_Control_Isolation.csv` 和 `PHASE3_2_V1_MESH_FAMILY_FAILURE_CN.md`。这些是 mesh-family fail-before，不是 radiation/diffraction/ABC/traction/A-B/RAO implementation failure；v1 没有运行 BEM。

## 2. v2 协同非嵌套 mesh family

v2 只改变离散控制，不改变任何物理量、outer extent 规则、production formula 或 v1 threshold。各级 theta count 精确等于该 body division 所产生的 legacy outer ring count，消除不必要的 count-transition sliver triangles。

```ini
NESTED_MESH = NO
MESH_FAMILY = coordinated_non_nested_L1_L2_L3_v2
```

| 控制量 | MESH_L1 | MESH_L2 | MESH_L3 |
|---|---:|---:|---:|
| body cube-face division | 4 | 6 | 8 |
| FS radial `[near,sponge]` | `[4,2]` | `[6,3]` | `[8,4]` |
| bottom radial `[near,sponge]` | `[4,2]` | `[6,3]` | `[8,4]` |
| configured wall vertical count | 3 | 4 | 6 |
| common FARFIELD theta count | 32 | 48 | 64 |
| FS outer theta count | 32 | 48 | 64 |
| bottom outer theta count | 32 | 48 | 64 |
| mesh transition radius | 20 m | 20 m | 20 m |
| sponge start radius | 20 m | 20 m | 20 m |

联合标称尺寸仍采用 v1 定义的七个倒数控制量几何平均：

```text
eta_L1 = 0.0951564092484917
eta_L2 = 0.0645140474965508
eta_L3 = 0.0475782046242459
r21 = 1.47497193155613
r32 = 1.35595800652962
```

## 3. 未改变的 gate

v1 的以下内容逐字按原语义继承：固定物理/1.5R outer extent；component mean/median h 与 transition/global h 记录；geometry、edge、Euler/manifold、normal、waterline、self-intersection 和 quality threshold；runtime provenance；raw/scaled rcond 与 residual gate；A33/B33/F3/RAO3/D3 scale floor；complex comparison；adjacent sensitivity；generalized unequal-ratio Richardson/GCI；`1.30:0.02:1.50 rad/s` resonance tracking；资源降级；production recommendation；Task 7 前置条件；所有禁止的事后调整。

```ini
THRESHOLD_RELAXATION_FROM_V1 = NO
PRODUCTION_PHYSICS_CHANGE = NO
V1_FAIL_PRESERVED = YES
V2_ALL_AFFECTED_CASES_MUST_BE_RERUN = YES
```
