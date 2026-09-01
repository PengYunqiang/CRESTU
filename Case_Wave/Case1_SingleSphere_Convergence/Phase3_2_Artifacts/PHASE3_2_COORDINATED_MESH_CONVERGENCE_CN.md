# CRESTU-1F Phase 3.2 协同网格收敛报告

```ini
PHASE3_2_OVERALL = BLOCKED
RUNTIME_PROVENANCE_GATE = PASS
Q_ONLY_LEGACY_GEOMETRY_REGRESSION = PASS_EXACT_HASH
MESH_CONTROL_ISOLATION = FAIL
MESH_GEOMETRY_GATE = FAIL
ALGEBRAIC_GATE = NOT_RUN_BLOCKED
SENTINEL_CONVERGENCE = NOT_RUN_BLOCKED
FORMAL_CONVERGENCE = NOT_RUN_BLOCKED
RESONANCE_CONVERGENCE = NOT_RUN_BLOCKED
ASYMPTOTIC_RANGE = NOT_DEMONSTRATED
PRODUCTION_MESH_RECOMMENDED = NO
READY_FOR_QUANTITATIVE_WAMIT_AUDIT = NO
PHYSICAL_VALIDATION = PARTIAL
READY_FOR_PRODUCTION = NO
```

## 已证明

MATLAB R2023b runtime/function provenance 通过；历史 q-only 1.5R geometry 在 5 个频率的 component/combined hash 完全一致。v4 生成 18 个 mesh-only cases，panel count 范围 2432–8792。同频率 outer extent 固定，原仓库路径未参与运行。

## Blocker

v1–v4 的原始 FAIL 均保留。v4 有 36 个 component quality rows 未通过；低频 FS 最大 aspect ratio 达 122.202，bottom/FS 最小角低于冻结 35 deg。0.5 rad/s transition median h 也未形成严格三层下降。根因已落到 base polar FS/bottom topology，继续处理需要另行授权的 generator redesign。

v1–v4 的原始 FAIL 均保留。v4 有 4.778345e+00 个 component quality rows 未通过；低频 FS 最大 aspect ratio 达 ## 未运行项目

按 fail-early 规则，没有组装 Rankine matrix，没有产生 A33/B33/F3/RAO3/D3、condition/residual、resonance 或 formal response。对应 CSV 明确写 `NOT_RUN_BLOCKED`；没有 figures，因为不存在可合法绘制的新响应数据。

## WAMIT

历史外部 IRR3 reference path 已恢复并 hash，原始数据再分发权限未确认。由于 mesh preconditions 失败，convention/discrepancy Task 7 未开始；WAMIT quantitative attribution 为 BLOCKED。
