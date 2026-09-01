# Phase 3.2A Topology Design Decision v2（冻结）

冻结时间：2026-09-01；在实现 connector-gap v2 和查看其 endpoint/v5 结果之前。

原始 `Phase3_2A_Topology_Design_Decision.md` 永久保留。initial v5 在 omega=1.0 停止，唯一失败子 gate 是 `all_component_mean_median_h_strictly_decrease`：bottom median h 为 `5.902969 → 6.476918 → 4.580403 m`。quality、mean h、counts、transition h、topology、outer extent、source/design hash 和 no-stale-reuse 全部 PASS。

Panel-population audit 证明 L2 median 落在固定 1:2 connector 的 minimum h。structured 5-block core median 为 `0.481 m` 且不是 quality blocker。

## 唯一设计变更

```text
v1: dr_1to2 = min(0.5*legacy_ring_chord, realized_outer_radial_step)
v2: dr_1to2 = min(0.4*legacy_ring_chord, realized_outer_radial_step)
```

其他 predeclared target、integer counts、radial layer schedule、hard quality gate 和 physics 全部不变。解析三角模板 oracle 在 normalized radius 50 m 给出：

| N inner | transition gap [m] | min angle [deg] | max AR |
|---:|---:|---:|---:|
| 32 | 3.9206856132 | 37.5229803271 | 1.6398439365 |
| 40 | 3.1383638291 | 37.7564036813 | 1.6319121662 |
| 56 | 2.2428178895 | 38.0195611902 | 1.6229212093 |

全部仍严格通过 35 deg 和 AR 12 gate。必须重新运行：

1. legacy q-only 5/5 exact regression；
2. omega=0.5/2.0 × L1/L2/L3 endpoint prototypes；
3. endpoint PASS 后完整 18-case v5 family。

```ini
THRESHOLD_RELAXATION = NO
COUNT_SCHEDULE_CHANGE = NO
PHYSICS_CHANGE = NO
LEGACY_DEFAULT_CHANGE = NO
INITIAL_V5_FAIL_PRESERVED = YES
```
