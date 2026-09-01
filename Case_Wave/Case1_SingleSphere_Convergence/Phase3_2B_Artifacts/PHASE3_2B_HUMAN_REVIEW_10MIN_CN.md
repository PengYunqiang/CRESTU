# Phase 3.2B 十分钟人工复核

1. 查看 `Phase3_2B_Final_Status.csv`，确认 overall=FAIL、production=NO。
2. 查看 frozen manifest：18 行、V2 显式、baseline reproducibility 全 PASS。
3. 查看 sentinel response：18 unique case、18 cache MISS、18 algebraic PASS。
4. 查看 condition/residual：最小 rcond 仍高于 1e-12，最大 residual 低于 1e-10。
5. 查看 convergence：24 mandatory 中 16 FAIL，floor sensitivity 均稳定。
6. 对照两张 sentinel 图，确认大差异未在报告中隐藏。
7. 查看 resonance/formal CSV，只能是 NOT_RUN_BLOCKED 且无数值结果。
8. 查看 legacy regression：START/END 各 8/8 PASS。
9. 查看 independent review：blocker 数与 final status 一致。
10. 查看 Git diff/hygiene：无 MAT/BMF/cache/temp/WAMIT raw，production source 未修改。
