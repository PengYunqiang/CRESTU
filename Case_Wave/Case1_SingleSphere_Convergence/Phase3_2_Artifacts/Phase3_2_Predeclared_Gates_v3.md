# CRESTU-1F Phase 3.2 预声明网格与判据（v3，冻结）

冻结时间：2026-09-01（在查看任何 v3 网格结果之前）

## 1. 继承与版本边界

v1/v2 gate、v1 全部 FAIL CSV 和 v2 fail-before report 永久保留。v3 不放宽 v1 的任何 geometry、quality、algebraic、adjacent sensitivity、resonance、Richardson/GCI、resource 或 production-readiness threshold；固定物理和按频率 1.5R outer extent 规则不变。

## 2. v3 mesh controls

| 控制量 | MESH_L1 | MESH_L2 | MESH_L3 |
|---|---:|---:|---:|
| body cube-face division | 4 | 6 | 8 |
| FS radial `[near,sponge]` | `[4,2]` | `[6,3]` | `[8,4]` |
| bottom radial `[near,sponge]` | `[4,2]` | `[6,3]` | `[8,4]` |
| configured wall vertical count | 3 | 4 | 6 |
| final FS/bottom/FARFIELD theta count | 64 | 96 | 128 |
| legacy ring count | 32 | 48 | 64 |
| mesh/sponge transition radius | 20 m | 20 m | 20 m |

```text
eta_L1 = 0.0707009403590185
eta_L2 = 0.0479337530744922
eta_L3 = 0.0353504701795093
r21 = 1.47497193155613
r32 = 1.35595800652962
```

## 3. Explicit-theta 渐进 transition policy

只在 PARA13 明确提供 theta controls 且 target count 大于 legacy count 时启用：

```text
minimum_transition_layers = ceil(log(N_final/N_legacy)/log(1.25))
radial_layer_count = max(ceil(extension/radial_target), minimum_transition_layers)
N_i = multiple-of-4 rounded linear interpolation from N_legacy to N_final
N_i/N_(i-1) <= 1.25
N_last = N_final
```

所有 count transition panel 仍属于 horizontal FS/bottom annulus，不进入 ABC-active FARFIELD wall。Wall 继续使用一个 common equal-count structured layout；outer ABC、normal/operator 公式不变。

q-only、legacy-absent 和未配置 theta control 的路径不得改变任何 vertices/count/hash。必须用冻结 Phase 3.1 representative geometry 做 exact hash regression。若 regression 不通过，v3 立即 FAIL，不运行 sentinel BEM。

## 4. v3 必须执行的证据链

1. v1 mesh-quality fail-before 与 v2 generator fail-before 保留；
2. production patch 仅限 `apply_rankine_outer_truncation.m` explicit-theta transition；
3. MATLAB Code Analyzer；
4. q-only/legacy geometry exact-hash regression；
5. v3 全部 18 个 mesh-only case 重跑并通过原 geometry threshold；
6. 只有 1–5 通过后才运行 clean sentinel BEM；
7. reviewer 独立核对 production diff 和 v1/v2 原始 FAIL。

```ini
THRESHOLD_RELAXATION_FROM_V1 = NO
PRODUCTION_PHYSICS_CHANGE = NO
PRODUCTION_MESH_EXPLICIT_CONTROL_PATCH = YES
LEGACY_Q_ONLY_GEOMETRY_MUST_BE_EXACT = YES
```
