function result = solve_uniform_flow(mesh_body, U_inf)
    np      = mesh_body.n_panels;
    centers = mesh_body.centers;
    normals = mesh_body.normals;
    verts   = mesh_body.vertices;
    e1_all  = mesh_body.e1;
    e2_all  = mesh_body.e2;

    fprintf('>>> 正在组装 BEM 影响矩阵 (N = %d)...\n', np);
    tic;

    A = zeros(np, np);
    b = -normals * U_inf(:); % 向量化右端项

    % 缓存 3D 诱导速度影响矩阵，避免后处理二次遍历
    Vx_mat = zeros(np, np);
    Vy_mat = zeros(np, np);
    Vz_mat = zeros(np, np);

    parfor i = 1:np
        x_i = centers(i, :);
        n_i = normals(i, :);
        
        row_A  = zeros(1, np);
        row_Vx = zeros(1, np);
        row_Vy = zeros(1, np);
        row_Vz = zeros(1, np);
        
        row_A(i) = 0.5; % 自感应跳跃项

        for j = 1:np
            if i ~= j
                P_j = squeeze(verts(j, :, :));
                [u_ij, v_ij, w_ij] = hess_smith_panel_velocity(P_j, x_i, ...
                                        centers(j, :), e1_all(j, :), e2_all(j, :), normals(j, :));
                
                row_Vx(j) = u_ij;
                row_Vy(j) = v_ij;
                row_Vz(j) = w_ij;
                row_A(j)  = u_ij * n_i(1) + v_ij * n_i(2) + w_ij * n_i(3);
            end
        end
        
        A(i, :)      = row_A;
        Vx_mat(i, :) = row_Vx;
        Vy_mat(i, :) = row_Vy;
        Vz_mat(i, :) = row_Vz;
    end
    fprintf('>>> 影响矩阵组装完成, 耗时: %.2f 秒\n', toc);

    % 线性求解源强
    sigma = A \ b;

    % 后处理：通过矩阵乘法直接得到诱导速度 (耗时 < 0.01s)
    V_ind = [Vx_mat * sigma, Vy_mat * sigma, Vz_mat * sigma];
    V_total = repmat(U_inf, np, 1) + V_ind + 0.5 * (sigma .* normals);

    U_mag2 = dot(U_inf, U_inf);
    Cp = 1.0 - sum(V_total.^2, 2) / U_mag2;

    result.sigma   = sigma;
    result.V_total = V_total;
    result.Cp      = Cp;
end