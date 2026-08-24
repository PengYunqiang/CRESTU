clear; clc; close all;

% =========================================================================
% 1. 读取主控制参数
% =========================================================================
cfg = read_config('CRESTU.cfg');

% =========================================================================
% 2. 组装单体 / 多体物面系统并变换至全局坐标系
% =========================================================================
body_list = cell(cfg.n_bodies, 1);
total_body_panels = 0;
for b = 1:cfg.n_bodies
    b_cfg = cfg.bodies(b);
    generate_body_bmf(b_cfg.mesh_file, 10.0, cfg.isx, cfg.isy, 5);
    mesh_loc = read_bmf(b_cfg.mesh_file);
    mesh_glob = transform_body_mesh(mesh_loc, b_cfg);
    body_list{b} = mesh_glob;
    total_body_panels = total_body_panels + mesh_glob.n_panels;
end

% 3. 提取全局水线 (Waterline)
waterline = extract_waterline(body_list{1}, cfg.z_tol);

% 4. 自动生成自由水面并导出 BMF (Type = 2)
mesh_fs = generate_free_surface_bmf(cfg.files.fs, waterline, cfg);

% 5. 处理有限水深分支 (海底 Type=4, 远场柱面 Type=5)
has_seabed = (cfg.water_depth > 0);
mesh_seabed = [];
mesh_farfield = [];
if has_seabed
    % 生成完整封闭海底网格
    mesh_seabed = generate_seabed_mesh(waterline, cfg, mesh_fs);
    write_bmf(cfg.files.seabed, mesh_seabed);
    mesh_farfield = generate_farfield_mesh(waterline, cfg, cfg.fs.nz_farfield);
    write_bmf(cfg.files.farfield, mesh_farfield);
    total_dofs = total_body_panels + mesh_fs.n_panels + ...
        mesh_seabed.n_panels + mesh_farfield.n_panels;
else
    total_dofs = total_body_panels + mesh_fs.n_panels;
end

% 6. 统计与封装
stats.total_body_panels = total_body_panels;
stats.fs_panels         = mesh_fs.n_panels;
stats.seabed_panels     = ternary(has_seabed, mesh_seabed.n_panels, 0);
stats.farfield_panels   = ternary(has_seabed, mesh_farfield.n_panels, 0);
stats.total_dofs        = total_dofs;

domain.cfg       = cfg;
domain.body_list = body_list;
domain.fs        = mesh_fs;
domain.seabed    = mesh_seabed;
domain.farfield  = mesh_farfield;
domain.waterline = waterline;
domain.stats     = stats;

% 打印流体域装配报告
fprintf('================= 全域 BMF 边界装配统计 =================\n');
fprintf(' 算例名称        : %s\n', cfg.case_name);
fprintf(' 浮体数量 (NBODY): %d (总物面面元: %d)\n', cfg.n_bodies, total_body_panels);
fprintf(' 自由水面面元数  : %d (截断外径 R_out = %.1f m)\n', mesh_fs.n_panels, cfg.fs.r_outer);
if has_seabed
    fprintf(' 海底边界面元数  : %d (水深 h = %.1f m)\n', mesh_seabed.n_panels, cfg.water_depth);
    fprintf(' 远场侧壁面元数  : %d\n', mesh_farfield.n_panels);
end
fprintf(' --------------------------------------------------------\n');
fprintf(' 单象限总未知数 (Total System DOFs): %d\n', total_dofs);
if ~isempty(body_list{1}.hydrostatics)
    fprintf(' 排水体积 Vz: %.4f m^3, 浮心 Zb: %.4f m\n', ...
        body_list{1}.hydrostatics.Vz, body_list{1}.hydrostatics.center_of_buoyancy(3));
end
fprintf('========================================================\n\n');

% 一键流体域三维可视化
plot_bmf_domain(domain,'wireframe');

% =========================================================================
% 7. 提取全域几何数组与物面专属数组
% =========================================================================
centers = []; normals = []; verts = []; body_areas = [];
for b = 1:cfg.n_bodies
    centers = [centers; domain.body_list{b}.centers];
    normals = [normals; domain.body_list{b}.normals];
    verts   = [verts; domain.body_list{b}.vertices];
    body_areas = [body_areas; domain.body_list{b}.areas];
end
centers = [centers; domain.fs.centers];
normals = [normals; domain.fs.normals];
verts   = [verts; domain.fs.vertices];
if has_seabed
    centers = [centers; domain.seabed.centers];
    normals = [normals; domain.seabed.normals];
    verts   = [verts; domain.seabed.vertices];
    
    centers = [centers; domain.farfield.centers];
    normals = [normals; domain.farfield.normals];
    verts   = [verts; domain.farfield.vertices];
end

body_centers = centers(1:domain.stats.total_body_panels, :);
body_normals = normals(1:domain.stats.total_body_panels, :);
nj = compute_generalized_normals(body_centers, body_normals, domain.body_list);

% =========================================================================
% 8. 频域求解与载荷后处理主循环
% =========================================================================
n_freqs = length(cfg.freq.omegas);
wave_headings = cfg.wave.headings; % 从 cfg 中安全读取浪向角数组
n_headings = length(wave_headings);

results(n_freqs) = struct();
pot_file_cache = sprintf('%s_Potential_Cache.mat', cfg.case_name);

for i_f = 1:n_freqs
    omega = cfg.freq.omegas(i_f);
    fprintf('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n');
    fprintf(' 开始处理频率工况 %d/%d: omega = %.4f rad/s\n', i_f, n_freqs, omega);
    fprintf('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n');
    
    phi_rad_body = [];
    phi_diff_body = [];
    phi_I_body = [];
    
    % --- 8.1 速度势求解阶段 (受 IPOTEN 控制) ---
    if cfg.run.ipoten == 1
        [A, S_body] = assemble_rankine_matrix(total_dofs, centers, normals, verts, domain.stats, omega, cfg.grav, cfg);
        
        fprintf('>>> 正在进行系统矩阵 LU 分解 (O(N^3))...\n');
        tic; 
        [L, U, P] = lu(A); 
        fprintf('>>> LU 分解完成, 耗时: %.2f 秒\n', toc);

        % 辐射势求解
        if cfg.run.irad == 1
            fprintf('>>> 极速求解辐射势 (全模态)...\n');
            b_rad = S_body * (1i * omega * nj);
            phi_rad = U \ (L \ (P * b_rad));
            phi_rad_body = phi_rad(1:domain.stats.total_body_panels, :);
        end

        % 绕射势求解
        if cfg.run.idiff == 1
            fprintf('>>> 极速求解绕射势 (全浪向)...\n');
            b_diff = complex(zeros(total_dofs, n_headings));
            phi_I_body = complex(zeros(domain.stats.total_body_panels, n_headings));
            for w = 1:n_headings
                [phi_I, dphi_I_dn] = compute_incident_wave(body_centers, body_normals, ...
                                        omega, cfg.grav, cfg.water_depth, wave_headings(w), 1.0);
                phi_I_body(:, w) = phi_I;
                b_diff(:, w) = S_body * (-dphi_I_dn);
            end
            phi_diff = U \ (L \ (P * b_diff));
            phi_diff_body = phi_diff(1:domain.stats.total_body_panels, :);
        end
        
        % 保存当前频率势函数到缓存 (可选按频点保存或全量保存)
        % 实际工程中可针对当前 frequency 保存至本地
    else
        fprintf('>>> [IPOTEN = 0] 跳过速度势求解，正在尝试从本地缓存加载...\n');
        if exist(pot_file_cache, 'file')
            load(pot_file_cache, 'saved_pots');
            phi_rad_body = saved_pots(i_f).phi_rad_body;
            phi_diff_body = saved_pots(i_f).phi_diff_body;
            phi_I_body = saved_pots(i_f).phi_I_body;
        else
            error('未找到速度势缓存文件 %s，请先设置 IPOTEN = 1 运行一次！', pot_file_cache);
        end
    end
    
    % --- 8.2 载荷与附加质量/阻尼后处理阶段 (受 IFORCE 控制) ---
    A_mat = []; B_mat = []; F_exc = [];
    if cfg.run.iforce == 1
        fprintf('>>> 正在进行水动力载荷与系数后处理...\n');
        
        % 辐射：附加质量与阻尼
        if cfg.run.irad == 1 && ~isempty(phi_rad_body)
            [A_mat, B_mat] = compute_hydrodynamic_coeffs(phi_rad_body, nj, body_centers, body_areas, omega, cfg.grav, cfg.rho);
        end
        
        % 绕射：激振力
        if cfg.run.idiff == 1 && ~isempty(phi_diff_body)
            F_exc = complex(zeros(6 * cfg.n_bodies, n_headings));
            for w = 1:n_headings
                if isempty(phi_I_body)
                    [phi_I, ~] = compute_incident_wave(body_centers, body_normals, ...
                                    omega, cfg.grav, cfg.water_depth, wave_headings(w), 1.0);
                    phi_I_col = phi_I;
                else
                    phi_I_col = phi_I_body(:, w);
                end
                
                phi_total = phi_I_col + phi_diff_body(:, w);
                pressure = 1i * omega * cfg.rho * phi_total;
                for dof = 1:(6 * cfg.n_bodies)
                    F_exc(dof, w) = sum(pressure .* nj(:, dof) .* body_areas);
                end
            end
        end
    else
        fprintf('>>> [IFORCE = 0] 跳过载荷计算。\n');
    end
    
    % 8.3 封装结果
    results(i_f).omega = omega;
    results(i_f).A_mat = A_mat;
    results(i_f).B_mat = B_mat;
    results(i_f).F_exc = F_exc;
    
    fprintf('>>> 频率 %.4f rad/s 处理完毕。\n', omega);
end

fprintf('\n================== 频域主程序执行完毕 ==================\n');

% 工具函数
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end