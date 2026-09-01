# CRESTU Phase 3.2A Base Polar Topology Remediation

```ini
PHASE3_2A_OVERALL = PASS_WITH_LIMITATIONS
QUALITY_METRIC_UNIT_TESTS = PASS_6_OF_6
V4_QUALITY_RECONCILIATION = PASS_72_OF_72
LEGACY_Q_ONLY_REGRESSION = PASS_5_OF_5_EXACT
ENDPOINT_COMPONENT_QUALITY = PASS_24_OF_24
V5_COMPONENT_QUALITY = PASS_72_OF_72
V5_CONTROL_ISOLATION = PASS_6_OF_6_FREQUENCIES
BEM_BRIDGE = PASS_2_OF_2
REVIEWER_BLOCKERS = 0
READY_FOR_PHASE3_2B_SENTINEL = YES
PHYSICAL_VALIDATION = PARTIAL
READY_FOR_PRODUCTION = NO
```

## 已修复内容

新增默认关闭的 `QUALITY_CONTROLLED_V2`。FS 使用无 smoothing annular O-grid；outer transition 使用固定 1:2 三角模板后接 equal-count quads；bottom 使用相同 annular rings 和已有 finite 5-block core，side radial layers 明确为 4/7/12。未引入 single-center fan。

## Fail-before 与版本证据

metric fixtures 排除 evaluator bug；FS AR=122.20 定位到 legacy base radial grading/relaxation，bottom 18–19 deg 定位到 zipper triangles。initial v5 与 gap-v2 median failures 均保留；最终 design-v3 在不改变 boundary nodes 的情况下修复 bottom median refinement。

## Mesh-only 验证

Endpoint 为 24 rows、完整 v5 为 18 manifest/72 quality rows；all quality、requested/realized h、strict counts、outer extent、manifold、orientation、waterline 和 self-intersection gates PASS。最大 v5 unknown count（仅成本 proxy）为 9744。

## 有限 BEM bridge

LEGACY N=2752，response 最大相对回归差 4.28e-13；新 mode N=4480，raw/scaled rcond=2.4e-05/5.87e-05，最大 algebraic/linear residual 3.06e-14。两者 cache MISS，A33/B33/F3/RAO3/D3 全有限。

## 限制

本阶段没有运行 formal sweep、resonance scan、WAMIT quantitative attribution，也没有证明 mesh convergence 或 physical production validation。新/旧响应差异不是 WAMIT 拟合结论。`READY_FOR_PHASE3_2B_SENTINEL=YES` 只授权下一阶段 sentinel。
