# CRESTU Phase 3.2A Polar Topology Design Decision（冻结）

冻结时间：2026-09-01；在实现 `QUALITY_CONTROLLED_V2`、查看 endpoint prototype 结果之前。

## 1. Fail-before 结论

独立 analytic fixtures 6/6 PASS；新 evaluator 从 mesh 逐行复算 blocked-v4 的 72/72 quality rows，得到相同的 36 PASS / 36 FAIL。因此 `METRIC_BUG=NO`。

- FS `AR=122.201862` 位于 MESH_L2、omega=0.5、base polar ring 8、sector 40、panel 320。该面板内侧两个顶点被 radial grading/legacy elliptic relaxation 拉回约 `r=20 m`，外侧两个顶点在 `r=350 m`。这是 base FS radial scheduler/relaxation 的真实几何缺陷。
- bottom 的 18–19 deg minimum angle 位于 outer-truncation 的 zipper triangles，不是 single-center fan。当前 bottom producer 已使用 adaptive 5-block structured core。
- omega=0.5 的 combined transition median h 为 L1/L2/L3=`35.6223/12.2430/15.5471 m`。L2 在 `r=20 m` 附近产生异常小面元，改变合并样本的中位数排序；不是 CSV 统计错误。

这些是多个已确认机制，不合并伪装成一个唯一 root cause。总体仍标记：

```ini
BASE_POLAR_TOPOLOGY_CLASSIFICATION = PRIMARY_DESIGN_LIMITATION
UNIQUE_ROOT_CAUSE = NO
```

## 2. Mode contract

```ini
POLAR_TOPOLOGY_MODE = LEGACY
POLAR_TOPOLOGY_MODE = QUALITY_CONTROLLED_V2
DEFAULT_MODE = LEGACY
```

未显式提供 `POLAR_TOPOLOGY_MODE` 的 production/config/历史 runner 必须继续使用 `LEGACY`，并保持 5 个冻结 q-only 1.5R component/combined hashes 完全一致。新模式只在 PARA13 显式请求时启用。

## 3. v5 coordinated targets

解析半球半径 `R=5 m`，mesh transition radius `Rt=20 m`。各级沿用 body divisions 4/5/7，并由 body waterline 精确给出 `N=32/40/56` 个 inner nodes。FS inner ring 与 body waterline逐点完全一致。

| quantity | L1 | L2 | L3 |
|---|---:|---:|---:|
| body division | 4 | 5 | 7 |
| inner/base theta count N | 32 | 40 | 56 |
| final wall theta count 2N | 64 | 80 | 112 |
| near radial layers | 4 | 5 | 6 |
| configured minimum outer radial layers | 2 | 3 | 4 |
| wall vertical control | 3 | 4 | 5 |
| inner chord target [m] | 0.9801714033 | 0.7845909573 | 0.5607044724 |
| requested near radial step [m] | 3.75 | 3.00 | 2.50 |
| chord at Rt [m] | 3.9206856132 | 3.1383638291 | 2.2428178895 |
| maximum requested outer radial step [m] | 31.3654849055 | 25.1069106329 | 17.9425431159 |

定义：

```text
h_theta_inner = 2*R*sin(pi/N)
dr_near = (Rt-R)/N_near
h_theta_transition = 2*Rt*sin(pi/N)
dr_outer_requested_max = 8*h_theta_transition
N_outer_layers = max(N_outer_configured,
                     ceil((R_legacy-Rt)/dr_outer_requested_max))
```

requested refinement ratios：near `r21=1.25, r32=1.20`；inner/transition tangential `r21≈1.2492, r32≈1.3992`。所有 realized counts 必须严格 L1<L2<L3；所有 realized component mean/median h 必须严格 L1>L2>L3。requested/realized h、ratio 和 integer counts 必须写入 endpoint/v5 evidence。

## 4. FS QUALITY_CONTROLLED_V2

1. 使用解析 annular O-grid；不执行 smoothing/elliptic relaxation。
2. 从 body waterline 到 `Rt` 使用 N 个等角节点和预声明 near radial layers。
3. 从 `Rt` 到 frequency-local legacy radius 使用 N 个等角节点；outer radial layers 由上述 requested maximum step 计算。
4. base outer ring 保持 N nodes。outer truncation 只允许一个预定义 1:2 connector，把 N 变为 2N；禁止 arbitrary zipper。
5. 1:2 connector 对每个 inner sector 使用固定三角模板，transition gap：

```text
dr_1to2 = min(0.5*legacy_ring_chord, realized_outer_radial_step)
```

6. connector 后到 final 1.5R ring 全部为 equal-count 2N quads；outer ring 与 FARFIELD wall top 逐点相同。

## 5. Bottom QUALITY_CONTROLLED_V2

1. 不引入 single-center fan。
2. 保留已存在的 finite adaptive 5-block structured core，但必须独立通过 35 deg/AR/warp gate。
3. core 外使用由 QUALITY_CONTROLLED_V2 FS 复制并下移到 seabed 的 annular rings。
4. outer truncation 使用与 FS 相同的预定义 1:2 connector；最终 2N ring 与 wall bottom 逐点相同。

## 6. Hard gates（不变）

- minimum corner angle `>=35 deg`；body AR `<=4`；near AR `<=4.5`；outer/transition/wall AR `<=12`；skewness `<=55 deg`；warp `<=8 deg`；poor fraction `<=5%`。
- endpoint prototype 仅 omega=0.5/2.0 × L1/L2/L3，要求 24/24 component rows PASS、transition median h 严格下降、四区域 count 均变化、outer extent exact、manifold/orientation/intersection/waterline PASS。
- endpoint 未 PASS 不得运行其余 12 cases。
- v5 全量要求 72/72 component rows PASS；只有此后才允许 omega=1.5 的两个有限 BEM bridge cases。
- legacy q-only exact regression 必须始终 5/5 PASS。

## 7. Production modification boundary

允许修改：config mode parser/manifest、FS/bottom mesh generation 的显式新 mode、explicit new-mode 1:2 horizontal connector、cache/source fingerprint 和 Phase 3.2A harness。

禁止修改：radiation/diffraction RHS、traction、outer ABC operator、source orientation、A/B extraction、RAO equation、WAMIT reference、任何 quality threshold。`LEGACY` 数值路径必须保持 exact。
