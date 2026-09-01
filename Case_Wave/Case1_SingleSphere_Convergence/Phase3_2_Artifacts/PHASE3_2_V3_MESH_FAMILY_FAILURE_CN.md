# Phase 3.2 v3 mesh family 失败记录

```ini
MESH_FAMILY_VERSION = v3
Q_ONLY_LEGACY_GEOMETRY_REGRESSION = PASS
MESH_CONTROL_ISOLATION_GEOMETRY = PASS
MESH_GEOMETRY_GATE = FAIL
ALGEBRAIC_GATE = NOT_RUN
PRODUCTION_PHYSICS_CHANGED = NO
```

v3 的 explicit-theta-only 渐进 transition 没有改变 5 个冻结 Phase 3.1 q-only 1.5R mesh 的任何 component/combined hash。18 个 v3 mesh-only case 的四区域 count、h、hash、outer extent 和 physical/source isolation 均通过。

原始 geometry FAIL 有两类：

1. MESH_L2（body division=6）在全部频率固定出现 32 条 nonmanifold edges，Euler characteristic 为 -6；MESH_L1/L3 对应为 0 条和 Euler=2。该固定共变归类为 mesh-family inconsistency，v4 不再使用 division=6。
2. v3 用 4 个 transition layers 把 endpoint count 从 N 增至 2N，每层仍跳 8/12/16 nodes。最小角约 22–31 deg，未达到冻结 35 deg；部分 FS/bottom aspect/skewness 也失败。

所有 v3 原始行保存在 `Phase3_2_V3_*.csv`，未删除或放宽 threshold。
