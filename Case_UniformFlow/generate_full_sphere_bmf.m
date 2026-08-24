function mesh_body = generate_full_sphere_bmf(filename, diameter, isx, isy, N)
% GENERATE_FULL_SPHERE_BMF 生成纯四边形完整封闭球体物面网格 (Full Sphere) 并导出为 .bmf (Type = 1)
%
% 算法:
%   采用标准 6 面立方体等角正切球面投影 (Cubed Sphere with equiangular tangent projection)
%   全域无极点退化面元，各面元正交性与长宽比极佳。

    if nargin < 1 || isempty(filename), filename = 'full_sphere_body.bmf'; end
    if nargin < 2, diameter = 10.0; end
    if nargin < 3, isx = 0; end
    if nargin < 4, isy = 0; end
    if nargin < 5, N = 8; end

    R = diameter / 2.0;

    % 1. 构建完整立方体 6 个面 (坐标范围 [-1, 1] x [-1, 1])
    grid_lin = linspace(-1, 1, 2 * N + 1);
    [U, V] = meshgrid(grid_lin, grid_lin);
    ones_mat = ones(size(U));

    faces = cell(6, 1);
    faces{1} = { ones_mat,  U,         V        }; % +X 面
    faces{2} = {-ones_mat, -U,         V        }; % -X 面
    faces{3} = { U,         ones_mat,  V        }; % +Y 面
    faces{4} = {-U,        -ones_mat,  V        }; % -Y 面
    faces{5} = { U,         V,         ones_mat }; % +Z 顶面 (封闭上半球)
    faces{6} = {-U,         V,        -ones_mat }; % -Z 底面 (下半球)

    raw_panels = [];

    % 2. 遍历 6 个面并等角投影至球面
    for f = 1:6
        Xf = faces{f}{1};
        Yf = faces{f}{2};
        Zf = faces{f}{3};
        [nr, nc] = size(Xf);

        for r = 1:nr-1
            for c = 1:nc-1
                c1 = [Xf(r, c),     Yf(r, c),     Zf(r, c)];
                c2 = [Xf(r, c+1),   Yf(r, c+1),   Zf(r, c+1)];
                c3 = [Xf(r+1, c+1), Yf(r+1, c+1), Zf(r+1, c+1)];
                c4 = [Xf(r+1, c),   Yf(r+1, c),   Zf(r+1, c)];

                p1 = project_sphere(c1, R);
                p2 = project_sphere(c2, R);
                p3 = project_sphere(c3, R);
                p4 = project_sphere(c4, R);

                p_center = (p1 + p2 + p3 + p4) / 4.0;

                % 对称性截断过滤
                if isx == 1 && (p_center(1) < -1e-6), continue; end
                if isy == 1 && (p_center(2) < -1e-6), continue; end

                % 确保法向严格向外指向流体 (+z / 径向背离球心)
                v12 = p2 - p1;
                v14 = p4 - p1;
                if dot(cross(v12, v14), p_center) < 0
                    p_temp = p2; p2 = p4; p4 = p_temp;
                end

                raw_panels = cat(1, raw_panels, reshape([p1; p2; p3; p4], [1, 4, 3]));
            end
        end
    end

    % 3. 构建结构体并写入 BMF
    total_panels = size(raw_panels, 1);
    mesh_body.header     = sprintf('Pure-Quad Full Sphere D=%.2fm, ISX=%d, ISY=%d', diameter, isx, isy);
    mesh_body.ulen       = 1.0;
    mesh_body.panel_type = repmat(1, [total_panels, 1]); % 1 = BODY
    mesh_body.isx        = isx;
    mesh_body.isy        = isy;
    mesh_body.n_panels   = total_panels;
    mesh_body.vertices   = raw_panels;

    write_bmf(filename, mesh_body);
    mesh_body = read_bmf(filename);

    % 4. 几何与静水力自检
    V_exact = 4.0 / 3.0 * pi * R^3;
    if isx == 1, V_exact = V_exact / 2.0; end
    if isy == 1, V_exact = V_exact / 2.0; end

    fprintf('>>> 全球体网格生成完毕: 文件 [%s], 面元数: %d\n', filename, total_panels);
    fprintf('    理论排水体积: %.4f m^3, 数值积分体积: %.4f m^3 (误差: %.4f%%)\n', ...
            V_exact, mesh_body.hydrostatics.V_mean, ...
            abs(mesh_body.hydrostatics.V_mean - V_exact)/V_exact * 100);
end

% -------------------------------------------------------------------------
% 辅助函数: 等角正切投影 (保持面元接近正交均匀)
% -------------------------------------------------------------------------
function P = project_sphere(pt, R)
    x = tan(pt(1) * pi / 4);
    y = tan(pt(2) * pi / 4);
    z = tan(pt(3) * pi / 4);
    P = R * [x, y, z] / sqrt(x^2 + y^2 + z^2);
end