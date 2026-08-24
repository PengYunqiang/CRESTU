function report = check_mesh_quality(domain, criteria)
% CHECK_MESH_QUALITY Execute the documented check_mesh_quality operation.
%
% Syntax:
%   report = check_mesh_quality(domain, criteria)
%
% Inputs:
%   domain          : [struct] Assembled Rankine boundary domain and configuration.
%   criteria        : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   report          : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%

    if nargin < 2, criteria = struct(); end
    if ~isfield(criteria, 'body_max_ar'), criteria.body_max_ar = 4.0; end % Maximum body-panel aspect ratio.
    if ~isfield(criteria, 'fs_near_max_ar'), criteria.fs_near_max_ar = 4.5; end % Near-field FS aspect ratio.
    if ~isfield(criteria, 'fs_sponge_max_ar'), criteria.fs_sponge_max_ar = 12.0; end % Sponge aspect ratio.
    if ~isfield(criteria, 'max_skewness'), criteria.max_skewness = 55.0; end % Maximum orthogonality error, deg.
    if ~isfield(criteria, 'min_angle'), criteria.min_angle = 35.0; end % Minimum corner angle, deg.
    if ~isfield(criteria, 'max_angle'), criteria.max_angle = 145.0; end % Maximum corner angle, deg.
    if ~isfield(criteria, 'max_warp_deg'), criteria.max_warp_deg = 8.0; end % Maximum panel warping, deg.

    maximum_parts = domain.cfg.n_bodies + 3;
    parts = cell(maximum_parts, 1);
    names = cell(maximum_parts, 1);
    part_types = zeros(maximum_parts, 1);
    part_count = 0;
    % Type 1: Body, Type 2: FS, Type 4: Seabed, Type 5: Farfield
    for b = 1:domain.cfg.n_bodies
        part_count = part_count + 1;
        parts{part_count} = domain.body_list{b};
        names{part_count} = sprintf('Body_%d', b);
        part_types(part_count) = 1;
    end
    part_count = part_count + 1;
    parts{part_count} = domain.fs;
    names{part_count} = 'FreeSurface';
    part_types(part_count) = 2;
    if domain.cfg.water_depth > 0
        if ~isempty(domain.seabed)
            part_count = part_count + 1;
            parts{part_count} = domain.seabed;
            names{part_count} = 'Seabed';
            part_types(part_count) = 4;
        end
        if ~isempty(domain.farfield)
            part_count = part_count + 1;
            parts{part_count} = domain.farfield;
            names{part_count} = 'Farfield';
            part_types(part_count) = 5;
        end
    end
    parts = parts(1:part_count);
    names = names(1:part_count);
    part_types = part_types(1:part_count);

    r_inner = domain.cfg.fs.r_inner;

    fprintf('\n======================= Rankine BEM 网格质量诊断报告 =======================\n');
    fprintf('%-13s | %-6s | %-8s | %-8s | %-8s | %-8s | %-8s\n', ...
            '边界部件', '面元数', '最大长宽比', '最大斜角(°)', '最小角(°)', '最大翘曲(°)', '不良率');
    fprintf('----------------------------------------------------------------------------\n');

    total_panel_count = sum(cellfun(@(item) item.n_panels, parts));
    all_verts = zeros(total_panel_count, 4, 3);
    all_is_bad = false(total_panel_count, 1);
    panel_cursor = 0;

    for p = 1:length(parts)
        m = parts{p};
        np = m.n_panels;
        verts   = m.vertices; % [np x 4 x 3]
        centers = m.centers;
        ptype   = part_types(p);

        aspect_ratios = zeros(np, 1);
        skew_angles   = zeros(np, 1);
        min_angles    = zeros(np, 1);
        max_angles    = zeros(np, 1);
        warps         = zeros(np, 1);
        is_bad        = false(np, 1);

        for i = 1:np
            P = squeeze(verts(i, :, :)); % [4 x 3]
            
            v = [P(2,:)-P(1,:); P(3,:)-P(2,:); P(4,:)-P(3,:); P(1,:)-P(4,:)];
            lens = sqrt(sum(v.^2, 2));

            aspect_ratios(i) = max(lens) / max(min(lens), 1e-12);

            angles = zeros(4, 1);
            for k = 1:4
                k_prev = mod(k - 2 + 4, 4) + 1;
                e_in  = -v(k_prev, :) / max(lens(k_prev), 1e-12);
                e_out =  v(k, :)      / max(lens(k), 1e-12);
                angles(k) = acosd(max(-1.0, min(1.0, dot(e_in, e_out))));
            end
            min_angles(i) = min(angles);
            max_angles(i) = max(angles);

            skew_angles(i) = max(abs(angles - 90.0));

            n1 = cross(v(1,:), v(2,:)); n2 = cross(v(3,:), v(4,:));
            if norm(n1)>1e-12 && norm(n2)>1e-12
                warps(i) = acosd(max(-1.0, min(1.0, dot(n1, n2)/(norm(n1)*norm(n2)))));
            end

            r_c = norm(centers(i, 1:2));
            is_sponge_zone = (ptype == 2 || ptype == 4) && (r_c > r_inner);

            if is_sponge_zone
                ar_limit = criteria.fs_sponge_max_ar;
            elseif ptype == 1
                ar_limit = criteria.body_max_ar;
            else
                ar_limit = criteria.fs_near_max_ar;
            end

            if ptype == 5
                if aspect_ratios(i) > criteria.fs_sponge_max_ar || warps(i) > criteria.max_warp_deg
                    is_bad(i) = true;
                end
            else
                if aspect_ratios(i) > ar_limit             || ...
                   skew_angles(i)   > criteria.max_skewness || ...
                   min_angles(i)    < criteria.min_angle    || ...
                   max_angles(i)    > criteria.max_angle    || ...
                   warps(i)         > criteria.max_warp_deg
                    is_bad(i) = true;
                end
            end
        end

        bad_rate = sum(is_bad) / np * 100.0;
        fprintf('%-13s | %-6d | %-8.2f | %-8.2f | %-8.2f | %-8.2f | %-7.2f%%\n', ...
                names{p}, np, max(aspect_ratios), max(skew_angles), min(min_angles), max(warps), bad_rate);

        report.(names{p}).aspect_ratio = aspect_ratios;
        report.(names{p}).skewness     = skew_angles;
        report.(names{p}).warping      = warps;
        report.(names{p}).is_bad       = is_bad;
        report.(names{p}).bad_rate     = bad_rate;

        panel_rows = panel_cursor + (1:np);
        all_verts(panel_rows, :, :) = verts;
        all_is_bad(panel_rows) = is_bad;
        panel_cursor = panel_cursor + np;
    end
    fprintf('============================================================================\n\n');

    generate_bem_advice(report, criteria);

    plot_quality_diagnostic(all_verts, all_is_bad);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function generate_bem_advice(report, criteria)
% GENERATE_BEM_ADVICE Execute the documented generate_bem_advice operation.
%
% Syntax:
%   generate_bem_advice(report, criteria)
%
% Inputs:
%   report          : [documented value] Input required by the implemented function contract.
%   criteria        : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   None; the function performs the documented file, plot, or validation action.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    fprintf('>>> BEM 求解器前处理优化建议:\n');
    fnames = fieldnames(report);
    any_issue = false;
    for i = 1:length(fnames)
        r = report.(fnames{i});
        if r.bad_rate > 5.0
            any_issue = true;
            fprintf(' * [%s] 不良面元占比 %.1f%%。', fnames{i}, r.bad_rate);
            if max(r.aspect_ratio) > criteria.fs_sponge_max_ar
                fprintf(' 径向拉伸过陡，建议降低 sponge_ratio 或增加海绵层径向段数。\n');
            elseif max(r.skewness) > criteria.max_skewness
                fprintf(' 局部面元剪切较明显，可适当增加 Winslow/Laplace 平滑迭代步数。\n');
            else
                fprintf(' 存在较大翘曲或过小内角，请检查局部网格过渡。\n');
            end
        end
    end
    if ~any_issue
        fprintf(' * [全部通过] 网格质量优秀！物面近场正交性良好，海绵层过渡自然，完全满足 Rankine BEM 精度要求。\n\n');
    end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function plot_quality_diagnostic(verts, is_bad)
% PLOT_QUALITY_DIAGNOSTIC Execute the documented plot_quality_diagnostic operation.
%
% Syntax:
%   plot_quality_diagnostic(verts, is_bad)
%
% Inputs:
%   verts           : [documented value] Input required by the implemented function contract.
%   is_bad          : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   None; the function performs the documented file, plot, or validation action.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    figure('Color', 'w', 'Position', [120, 120, 1000, 720], ...
        'Name', 'BEM Mesh Quality Verification');
    ax = gca; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    view(ax, 135, 30);
    xlabel('X (m)', 'FontWeight', 'bold');
    ylabel('Y (m)', 'FontWeight', 'bold');
    zlabel('Z (m)', 'FontWeight', 'bold');
    title('Rankine BEM 网格质量诊断图 [绿色: 合格面元, 红色: 待优化面元]', 'FontSize', 12);

    X = squeeze(verts(:, :, 1))';
    Y = squeeze(verts(:, :, 2))';
    Z = squeeze(verts(:, :, 3))';

    good_idx = find(~is_bad);
    if ~isempty(good_idx)
        patch(ax, X(:, good_idx), Y(:, good_idx), Z(:, good_idx), 'w', ...
              'FaceColor', [0.88 0.94 0.88], 'EdgeColor', [0.45 0.55 0.45], ...
              'LineWidth', 0.2, 'FaceAlpha', 0.6);
    end

    bad_idx = find(is_bad);
    if ~isempty(bad_idx)
        patch(ax, X(:, bad_idx), Y(:, bad_idx), Z(:, bad_idx), 'r', ...
              'FaceColor', [0.95 0.25 0.25], 'EdgeColor', [0.6 0.1 0.1], ...
              'LineWidth', 1.2, 'FaceAlpha', 0.9);
    end
end
