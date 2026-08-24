# CRESTU 修改与新增文件清单

本清单记录本轮实际重构范围。由于工作目录不是 Git 仓库，无法依据提交历史区分原始文件与此前未跟踪改动；
以下 53 个模块文件均经过统一英文文档、注释和格式处理，其中标注“核心”的文件同时包含算法修改。

## 模块文件

- `1.Input/read_config.m`（核心：PARA 字符串、数组和多体质量属性解析）
- `2.Mesh/build_bmf_domain.m`（核心：多体形状安全合并和按频率环境入口）
- `2.Mesh/check_mesh_quality.m`
- `2.Mesh/complete_waterline_by_symmetry.m`
- `2.Mesh/expand_mesh_by_symmetry.m`
- `2.Mesh/extract_waterline.m`
- `2.Mesh/generate_body_bmf.m`
- `2.Mesh/generate_box_bmf.m`
- `2.Mesh/generate_farfield_mesh.m`
- `2.Mesh/generate_free_surface_bmf.m`
- `2.Mesh/generate_multibody_free_surface_bmf.m`
- `2.Mesh/generate_reduced_seabed_mesh.m`
- `2.Mesh/generate_seabed_disk_mesh.m`
- `2.Mesh/generate_seabed_mesh.m`
- `2.Mesh/merge_domain_geometry.m`（核心：列向量归一化与索引偏移）
- `2.Mesh/plot_bmf_domain.m`
- `2.Mesh/read_bmf.m`
- `2.Mesh/reduce_mesh_by_symmetry.m`
- `2.Mesh/transform_body_mesh.m`
- `2.Mesh/tune_sponge_layer.m`（核心：波长约束和动态 Rayleigh 参数）
- `2.Mesh/write_bmf.m`
- `3.HessSmith/hess_smith_panel_velocity.m`
- `4.Potential/assemble_rankine_matrix.m`（核心：parfor 切片安全、像源叠加）
- `4.Potential/compute_generalized_normals.m`
- `4.Potential/compute_incident_wave.m`
- `4.Potential/decompose_incident_wave_symmetry.m`
- `4.Potential/get_mode_parities.m`
- `4.Potential/load_potential_cache.m`
- `4.Potential/rankine_panel_integrals.m`
- `4.Potential/save_potential_cache.m`
- `4.Potential/solve_complex_system.m`
- `4.Potential/solve_radiation_freq.m`
- `4.Potential/solve_wave_dispersion.m`
- `4.Potential/symmetry_force_weights.m`
- `5.Force/assemble_mass_matrix.m`
- `5.Force/compare_wamit_results.m`
- `5.Force/compute_haskind_excitation.m`
- `5.Force/compute_hydrodynamic_coeffs.m`（核心：互易诊断和对称投影）
- `5.Force/compute_hydrostatic_matrix.m`
- `5.Force/compute_radiation_damping_energy.m`
- `5.Force/compute_wave_excitation.m`
- `5.Force/export_surface_pressure.m`
- `5.Force/plot_surface_pressure.m`
- `5.Force/read_wamit_excitation.m`
- `5.Force/read_wamit_first_order.m`
- `5.Force/read_wamit_rao.m`
- `5.Force/solve_rao.m`
- `6.MeanDriftLoads/compute_drift_farfield.m`（核心：Kochin 法向、配对和尺度）
- `6.MeanDriftLoads/compute_drift_nearfield.m`（核心：Pinkster 四项积分）
- `6.MeanDriftLoads/estimate_surface_kinematics.m`（核心：二次 MLS 梯度/Hessian）
- `6.MeanDriftLoads/expand_scalar_by_symmetry.m`
- `6.MeanDriftLoads/export_drift_loads.m`（核心：力/力矩无量纲尺度）
- `6.MeanDriftLoads/read_wamit_mean_drift.m`（核心：ULEN 与报告尺度）

## 主控、测试与公共脚本

- `Case_Wave/run_frequency_domain_case.m`（核心：频率局部外域、缓存版本 4、结果版本 5）
- `Case_Wave/Run_All_Validation.m`（新增）
- `Case_Wave/Case1_SingleSphere_Convergence/Run_SingleSphere_Convergence.m`（新增）
- `Case_Wave/Case2_TwoSpheres_Interaction/Run_TwoSpheres_Interaction.m`（新增）
- `Case_Wave/Common_Scripts/read_wamit_dataset.m`（新增）
- `Case_Wave/Common_Scripts/match_reference_grid.m`（新增）
- `Case_Wave/Common_Scripts/export_chinese_figure.m`（新增）
- `Case_Wave/Common_Scripts/Build_Validation_Tables.m`（新增）
- `Case_Wave/Common_Scripts/Plot_SingleBody_Results.m`（新增）
- `Case_Wave/Common_Scripts/Plot_MultiBody_Results.m`（新增）
- `Case_Wave/Common_Scripts/Plot_MeanDrift_Comparison.m`（新增）
- `Case_Wave/Benchmark_Symmetry.m`
- `Case_Wave/Run_Convergence_Study.m`
- `Case_Wave/Test_MultiBody_Wave.m`
- `Case_Wave/Validate_Hemisphere_WAMIT.m`

## 配置、网格和可重复输出

- `Case_Wave/` 下 22 个 `.cfg` 已统一为 16 频点和 `0/45/90 deg` 浪向。
- 新增 `Case1_SingleSphere_Convergence/Mesh_Coarse`、`Mesh_Medium`、`Mesh_Fine` 的独立配置、网格、缓存和结果。
- 新增 `Case2_TwoSpheres_Interaction` 的双体配置、两体网格、缓存和结果。
- 新增单体、双体对比 CSV/MAT，以及三幅 300 DPI 中文对比图。
- 新增 `tools/standardize_matlab_docs.ps1`、`tools/normalize_matlab_comments.ps1` 和
  `tools/update_case_frequency_grid.ps1`，用于重复执行文档和配置规范化。

未删除任何既有网格、转换器、缓存或测试脚本，`WAMIT/` 全程只读。
