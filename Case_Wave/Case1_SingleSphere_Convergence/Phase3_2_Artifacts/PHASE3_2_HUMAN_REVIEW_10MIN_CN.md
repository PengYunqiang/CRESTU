# Phase 3.2 十分钟人工复核

1. 先读 `Phase3_2_Final_Status.csv`，确认 overall BLOCKED、READY_FOR_PRODUCTION=NO。
2. 读 v1–v4 gate/report，确认 threshold 未放宽且旧 FAIL 未覆盖。
3. 在 `Phase3_2_Mesh_Manifest.csv` 核对 18 rows、四区域 count/hash、1.5R outer extent。
4. 在 `Phase3_2_Mesh_Quality.csv` 查看 36 个 FAIL rows，尤其 bottom center fan 和低频 FS grading。
5. 读 `Phase3_2_Legacy_Mesh_Regression.csv`，确认 q-only 5/5 exact PASS。
6. 检查四个 convergence/condition/resonance CSV 均明确 NOT_RUN_BLOCKED，而不是空白或伪造数值。
7. 读 production diff，确认只改 explicit-theta horizontal annulus transition，没有 ABC/RHS/traction/force/RAO patch。
8. 检查没有 BMF/cache/MAT/temp/WAMIT raw 新增文件。
