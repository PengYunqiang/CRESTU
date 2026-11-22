function summary = run_appendix_a(outputRoot)
%RUN_APPENDIX_A Reproduce the P2 Appendix A semicircle and rectangular checks.
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    if nargin < 1
        outputRoot = fullfile(root, 'hf_results');
    end
    assert(~isfolder(outputRoot) && ~isfile(outputRoot), ...
        'CRESTU2D:OutputExists', 'Choose a new output directory.');
    mkdir(outputRoot);
    [ok, outputFolder] = fileattrib(outputRoot);
    assert(ok);
    outputRoot = outputFolder.Name;
    previousFolder = pwd;
    restoreFolder = onCleanup(@() cd(previousFolder)); %#ok<NASGU>
    cd(root);
    CRESTU_SETUP('2d');
    gates = jsondecode(fileread(fullfile(root, 'research', 'p2_hf', 'ACCEPTANCE_CRITERIA.json')));

    opt = struct('rho', 1025, 'depth', Inf, 'reference', [0 0], 'quadrature_order', 12);
    panels = [32; 64; 128; 256];
    mass = zeros(size(panels));
    reciprocity = zeros(size(panels));
    for j = 1:numel(panels)
        theta = linspace(-pi, 0, panels(j) + 1)';
        nodes = [cos(theta), sin(theta)];
        nodes([1 end], 2) = 0;
        result = crestu2d.hf_solve(nodes, opt);
        mass(j) = result.A(2, 2);
        reciprocity(j) = result.reciprocity;
    end
    exact = repmat(opt.rho * pi / 2, size(panels));
    relative = abs(mass ./ exact - 1);
    circle = table(panels, mass, exact, relative, reciprocity, ...
        'VariableNames', {'panels', 'A', 'exact', 'relative', 'reciprocity'});
    writetable(circle, fullfile(outputRoot, 'analytic_circle.csv'));
    assert(relative(end) <= gates.analytic_semicircle_heave_relative_max, ...
        'CRESTU2D:AppendixCircle', 'Semicircle failed the analytical error criterion.');

    run_hf_validation('table', outputRoot);
    results = readtable(fullfile(outputRoot, 'table', 'results.csv'));
    coarse = sortrows(results(results.panels_per_edge == 64, :), {'width_m', 'draft_m'});
    fine = sortrows(results(results.panels_per_edge == 128, :), {'width_m', 'draft_m'});
    keys = {'width_m', 'draft_m'};
    assert(height(fine) == 38 && height(results) == 76 && ...
        isequal(coarse(:, keys), fine(:, keys)) && height(unique(fine(:, keys))) == 38, ...
        'CRESTU2D:AppendixGrid', 'Expected 38 paired rectangular nodes.');
    fine.BEM_N64_kg_per_m = coarse.A_heave_kg_per_m;
    fine.body_grid_change_relative = abs(fine.A_heave_kg_per_m ./ fine.BEM_N64_kg_per_m - 1);
    writetable(fine, fullfile(outputRoot, 'rectangle_refinement.csv'));
    assert(all(isfinite(fine.body_grid_change_relative)) && ...
        all(fine.body_grid_change_relative <= gates.body_finest_successive_heave_relative_max) && ...
        all(fine.linear_residual <= gates.linear_residual_max) && ...
        all(fine.reciprocity <= gates.normalized_matrix_reciprocity_max) && ...
        all(fine.minimum_scaled_eigenvalue >= gates.normalized_minimum_eigenvalue_min), ...
        'CRESTU2D:AppendixRectangle', 'At least one rectangular node failed a BEM check.');

    summary = struct('rectangle_nodes', height(fine), ...
        'semicircle_finest_A_heave_kg_per_m', mass(end), ...
        'semicircle_finest_relative_error', relative(end), ...
        'maximum_rectangle_grid_change_relative', max(fine.body_grid_change_relative));
    fprintf(['APPENDIX_A PASS | %d rectangle nodes | semicircle relative error %.10g | ', ...
        'maximum rectangle grid change %.10g\n'], summary.rectangle_nodes, ...
        summary.semicircle_finest_relative_error, summary.maximum_rectangle_grid_change_relative);
end
