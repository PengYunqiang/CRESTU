# Phase 3.2 v2 mesh family fail-before

```ini
MESH_FAMILY_VERSION = v2
FIRST_CASE = MESH_L1, omega=0.5 rad/s
MESH_GENERATION = FAIL_BEFORE
ERROR_IDENTIFIER = CRESTU:OuterTruncationThetaRounding
ALGEBRAIC_GATE = NOT_RUN
PRODUCTION_PHYSICS_CHANGED = NO
```

v2 试图让 endpoint theta count 与 legacy ring count 相同，以消除 v1 的单层 count-transition sliver triangles。Production generator 的冻结 resolution audit 正确地拒绝了这一设计：Phase 3.2 的 1.5R endpoint 半径约为 legacy radius 的两倍，保持原 count 会使 circumferential spacing 比 frozen chord target 粗超过 5%。该 fail 发生在首个 mesh build 内，未组装 Rankine matrix、未产生 cache、未查看响应。

结论：endpoint count 必须约随半径扩大为 2N，但 count transition 必须分布在多个非 ABC horizontal annulus layers 中，不能像 v1 一样单层从 N 跳到约 2N。v3 只为 explicit-theta 路径增加渐进 count transition；q-only/legacy 路径必须保持原几何/hash。
