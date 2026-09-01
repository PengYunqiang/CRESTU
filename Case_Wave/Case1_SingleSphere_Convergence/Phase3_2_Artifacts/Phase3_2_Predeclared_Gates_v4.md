# CRESTU-1F Phase 3.2 预声明网格与判据（v4，冻结）

冻结时间：2026-09-01（在查看任何 v4 网格结果之前）

## 1. 继承

v1/v2/v3 gate、FAIL CSV/report 和 v3 exact q-only legacy regression 永久保留。v4 继承 v1 的全部 threshold、fixed physics、1.5R outer extent、complex comparison、GCI、resonance、resource 和 production-readiness 规则，不放宽任何 gate。

## 2. v4 family

| 控制量 | MESH_L1 | MESH_L2 | MESH_L3 |
|---|---:|---:|---:|
| body cube-face division | 4 | 5 | 7 |
| FS radial `[near,sponge]` | `[4,2]` | `[5,3]` | `[6,4]` |
| bottom radial `[near,sponge]` | `[4,2]` | `[5,3]` | `[6,4]` |
| wall vertical count | 3 | 4 | 5 |
| legacy theta count | 32 | 40 | 56 |
| final FS/bottom/wall theta count | 64 | 80 | 112 |

```text
eta_L1 = 0.0707009403590185
eta_L2 = 0.0550177552570247
eta_L3 = 0.0412542094946932
r21 = 1.28505679718715
r32 = 1.33362766929522
```

## 3. v4 gradual transition

仅 explicit-theta 路径使用：

```text
minimum_transition_layers = (N_final-N_legacy)/4
N_i increases by at most 4 nodes per radial layer
radial_layer_count = max(ceil(extension/radial_target), minimum_transition_layers)
N_last = N_final
```

因此三层最少分别使用 8/10/14 个 non-ABC horizontal-annulus transition layers。q-only/legacy path 仍必须通过相同 exact hash regression。

```ini
THRESHOLD_RELAXATION = NO
V1_V2_V3_FAIL_PRESERVED = YES
LEGACY_Q_ONLY_GEOMETRY_MUST_BE_EXACT = YES
```
