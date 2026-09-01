function [A, S_body, assembly_info] = assemble_rankine_matrix(N, centers, normals, vertices, stats, omega, g, cfg, parity)
% ASSEMBLE_RANKINE_MATRIX Assemble the mixed-boundary Rankine boundary-integral system.
%
% Syntax:
%   [A, S_body, assembly_info] = assemble_rankine_matrix(N, centers, normals, vertices, stats, omega, g, cfg, parity)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   N                  - [scalar] Total number of boundary panels, dimensionless.
%   centers            - [N x 3] Panel collocation points in global coordinates, [m].
%   normals            - [N x 3] Unit panel normals, dimensionless.
%   vertices           - [N x 4 x 3] Ordered quadrilateral panel vertices, [m].
%   stats              - [struct] Boundary-component panel counts, dimensionless.
%   omega              - [scalar] Angular frequency, [rad/s].
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%   parity             - [1 x 2] Reflection parity signs for the x = 0 and y = 0 planes, dimensionless.
%
% Outputs:
%   A                  - [N x N] Complex boundary-integral influence matrix, dimensionless after boundary-condition scaling.
%   S_body             - [N x Nbpan] Body-source influence matrix, [m].
%   assembly_info      - [struct] Assembly timing, wavenumber, image, and sponge-layer diagnostics.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if nargin < 9 || isempty(parity)
        parity = [1, 1];
    end
    validateattributes(N, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(omega, {'numeric'}, {'scalar','real','positive','finite'});
    validateattributes(g, {'numeric'}, {'scalar','real','positive','finite'});
    parity = reshape(parity, 1, 2);
    if any(abs(parity) ~= 1)
        error('CRESTU:SymmetryParity','Parity entries must be +1 or -1.');
    end
    centers = reshape(centers, N, 3);
    normals = reshape(normals, N, 3);
    vertices = reshape(vertices, N, 4, 3);
    n_body = stats.total_body_panels;
    n_fs = stats.fs_panels;
    n_sb = stats.seabed_panels;
    n_ff = N - n_body - n_fs - n_sb;
    if n_ff < 0 || (isfield(stats,'farfield_panels') && stats.farfield_panels ~= n_ff)
        error('CRESTU:DomainCounts','Panel statistics do not match N=%d.', N);
    end
    if any(~isfinite(centers(:))) || any(~isfinite(normals(:))) || any(~isfinite(vertices(:)))
        error('CRESTU:NonfiniteGeometry','BEM geometry contains NaN or Inf.');
    end
    source_orientation = get_rankine_source_orientation( ...
        stats, normals, vertices);
    [image_vertices, image_normals, image_weights, image_labels] = prepare_images( ...
        vertices, normals, cfg.isx, cfg.isy, parity);

    [k0, ~] = solve_wave_dispersion(omega, g, cfg.water_depth);
    free_surface_bc = get_rankine_free_surface_robin_bc( ...
        centers, stats, omega, g, cfg);
    nu0 = free_surface_bc.nu0; % [1/m], exact finite-depth linear FS Robin scale
    free_surface_robin_convention = ...
        'exp(+i*omega*t):dphi/dz=(omega^2/g)*(1-i*mu(r))*phi';
    sponge_start = free_surface_bc.startRadius;
    sponge_width = free_surface_bc.width;
    effective_mu0 = free_surface_bc.actualMu0;
    sponge_manifest = struct('inputMode', 'legacy-direct-fallback');
    if isfield(cfg, 'phase2_2_controls') && ...
            isfield(cfg.phase2_2_controls, 'effective')
        sponge_state = cfg.phase2_2_controls.effective;
        sponge_start = sponge_state.spongeStartRadiusM;
        sponge_width = sponge_state.spongeWidthM;
        effective_mu0 = sponge_state.actualMu0;
        sponge_manifest = cfg.phase2_2_controls.manifest;
    end
    fs_idx = n_body + (1:n_fs);
    nu = zeros(N, 1);
% exp(i*w*t): exp(-i*k*r), hence Im(k)<0 gives spatial decay.
    nu(fs_idx) = free_surface_bc.nu;
    outer_bc = get_rankine_outer_absorbing_bc(centers, normals, stats, ...
        omega, g, cfg.water_depth, cfg, k0, effective_mu0);
    geometryHash = local_geometry_hash(centers, normals, vertices);
    if isfield(cfg, 'phase2_2_controls') && ...
            isfield(cfg.phase2_2_controls, 'effective') && ...
            isfield(cfg.phase2_2_controls.effective, 'mergedGeometrySHA256')
        geometryHash = cfg.phase2_2_controls.effective.mergedGeometrySHA256;
    end
    boundary_operator_state = get_rankine_boundary_operator_state( ...
        free_surface_bc.cacheSpecification, ...
        outer_bc.cacheSpecification, geometryHash);
    % <<<CORE>>> assemble_rankine_influence_matrix, paper_eq=green_third_identity, benchmark=single_sphere_grid_audit
    A = complex(zeros(N, N));
    S_body = complex(zeros(N, n_body));
    auditOuterColumn = n_ff > 0 && isfield(cfg, 'outer_truncation') && ...
        isfield(cfg.outer_truncation, 'effective') && ...
        cfg.outer_truncation.effective.geometryModified;
    auditOuterGlobalIndex = n_body + n_fs + n_sb + 1;
    auditOuterSingleLayer = complex(zeros(N, auditOuterColumn));
    auditOuterDoubleLayer = complex(zeros(N, auditOuterColumn));
    inv4pi = 1 / (4 * pi);
    fprintf('[INFO] Rankine assembly N=%d, parity=[%+d,%+d], images=%d\n', ...
        N, parity(1), parity(2), numel(image_weights)); timer = tic;
    for i = 1:N
        rowA = complex(zeros(1, N));
        rowS = complex(zeros(1, n_body));
        x = centers(i, :);
        for j = 1:N
            Gsum = 0;
            Dsum = 0;
            for im = 1:numel(image_weights)
                panel = reshape(image_vertices(j, :, :, im), 4, 3);
                normal = reshape(image_normals(j, :, im), 1, 3);
                [G, dGdn] = rankine_panel_integrals(panel, x, normal);
                Gsum = Gsum + image_weights(im) * G;
                Dsum = Dsum + image_weights(im) * dGdn;
            end
            Sij = Gsum * inv4pi;
% Only the double-layer source derivative is converted to the canonical
% inward-to-fluid source orientation. S, physical normals, right-hand sides,
% and the +0.5 self term retain their existing definitions.
            Dsum = source_orientation.columnSigns(j) * Dsum;
            Dij = -Dsum * inv4pi;
            if i == j
                Dij = Dij + 0.5;
            end
            if auditOuterColumn && j == auditOuterGlobalIndex
                auditOuterSingleLayer(i) = Sij;
                auditOuterDoubleLayer(i) = Dij;
            end
            if j <= n_body
                rowA(j) = Dij;
                rowS(j) = Sij;
            elseif j <= n_body + n_fs
                rowA(j) = Dij - nu(j) * Sij;
            elseif j <= n_body + n_fs + n_sb
                rowA(j) = Dij;
            else
% First-order local asymptotic absorbing BC with stored-inward normals.
% This source-local gamma applies to every collocation row; it is not an
% exact finite-radius Dirichlet-to-Neumann operator.
                outerSourceIndex = j - (n_body + n_fs + n_sb);
                rowA(j) = Dij + outer_bc.gamma(outerSourceIndex) * Sij;
            end
        end
        A(i, :) = rowA;
        S_body(i, :) = rowS;
    end
    elapsed = toc(timer);
    outerSourceColumnAudit = struct('enabled', false);
    if auditOuterColumn
        gammaAudit = outer_bc.gamma(1);
        reconstructedColumn = auditOuterDoubleLayer + ...
            gammaAudit * auditOuterSingleLayer;
        assembledColumn = A(:, auditOuterGlobalIndex);
        reconstructionResidual = norm(assembledColumn - reconstructedColumn) / ...
            max(norm(assembledColumn), eps);
        assert(reconstructionResidual <= 1e-14, ...
            'CRESTU:OuterTruncationColumnAudit', ...
            'Instrumented outer column is not D+gamma_j*S.');
        outerSourceColumnAudit = struct('enabled', true, ...
            'sourceGlobalIndex', auditOuterGlobalIndex, ...
            'sourceLocalIndex', 1, 'gamma', gammaAudit, ...
            'canonicalFormula', 'D+gamma_j*S', ...
            'reconstructionRelativeResidual', reconstructionResidual, ...
            'singleLayerColumn', auditOuterSingleLayer, ...
            'doubleLayerColumn', auditOuterDoubleLayer, ...
            'assembledColumn', assembledColumn, ...
            'singleLayerHash', hash_complex_column(auditOuterSingleLayer), ...
            'doubleLayerHash', hash_complex_column(auditOuterDoubleLayer), ...
            'assembledColumnHash', hash_complex_column(assembledColumn));
    end
    assembly_info = struct('elapsed_seconds', elapsed,'parity', parity, ...
'image_count', numel(image_weights),'image_labels', {image_labels}, ...
'wavenumber', k0,'wavelength', 2 * pi / k0,'effective_mu0', effective_mu0, ...
'sponge_start_radius', sponge_start, 'sponge_width', sponge_width, ...
'sponge_end_radius', sponge_start + sponge_width, ...
'sponge_control_manifest', sponge_manifest, 'n_unknowns', N, ...
'free_surface_robin_nu0', nu0, ...
'free_surface_robin_base_scale', 'omega^2/g', ...
'free_surface_robin_convention', free_surface_robin_convention, ...
'free_surface_robin_coefficients', free_surface_bc.cacheSpecification, ...
'boundary_operator_state', boundary_operator_state, ...
'outer_absorbing_boundary', outer_bc.cacheSpecification, ...
'double_layer_source_normal_convention', char(source_orientation.convention), ...
'double_layer_component_names', source_orientation.componentNames, ...
'double_layer_component_signs', source_orientation.componentSigns, ...
'double_layer_orientation_hash', source_orientation.signatureHash, ...
'double_layer_minimum_winding_alignment', ...
source_orientation.minimumWindingAlignment);
    assembly_info.outer_source_column_audit = outerSourceColumnAudit;
    if isfield(cfg, 'outer_truncation') && ...
            isfield(cfg.outer_truncation, 'effective')
        assembly_info.outer_truncation = cfg.outer_truncation.effective;
    else
        assembly_info.outer_truncation = struct('enabled', false, ...
            'geometryModified', false);
    end
    % <<</CORE>>>
    fprintf('[OK] Rankine assembly completed in %.3f s.\n', elapsed);
end

function hashText = local_geometry_hash(centers, normals, vertices)
% LOCAL_GEOMETRY_HASH Fallback identity for direct assembler calls.
    values = [centers(:); normals(:); vertices(:)];
    bytes = [reshape(typecast(uint64(size(values)), 'uint8'), 1, []), ...
        reshape(typecast(double(values), 'uint8'), 1, [])];
    hashText = sha256_hash(bytes);
end

function hashText = hash_complex_column(values)
% HASH_COMPLEX_COLUMN Hash exact real/imaginary column bytes for audit.
    bytes = [reshape(typecast(double(real(values(:))), 'uint8'), 1, []), ...
        reshape(typecast(double(imag(values(:))), 'uint8'), 1, [])];
    hashText = sha256_hash(bytes);
end

function [all_vertices, all_normals, weights, labels] = prepare_images(vertices, normals, isx, isy, parity)
% PREPARE_IMAGES Construct reflected panel images and parity weights for symmetry reduction.
%
% Syntax:
%   [all_vertices, all_normals, weights, labels] = prepare_images(vertices, normals, isx, isy, parity)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   vertices           - [N x 4 x 3] Ordered quadrilateral panel vertices, [m].
%   normals            - [N x 3] Unit panel normals, dimensionless.
%   isx                - [logical scalar] Symmetry-reduction flag for the x = 0 plane.
%   isy                - [logical scalar] Symmetry-reduction flag for the y = 0 plane.
%   parity             - [1 x 2] Reflection parity signs for the x = 0 and y = 0 planes, dimensionless.
%
% Outputs:
%   all_vertices       - [N x 4 x 3 x Nimage] Original and reflected panel vertices, [m].
%   all_normals        - [N x 3 x Nimage] Original and reflected unit panel normals, dimensionless.
%   weights            - [numeric array] Symmetry, interpolation, or quadrature weights, dimensionless.
%   labels             - [Nimage x 1 cell] Human-readable reflection-image labels.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    flags = [0, 0];
    if isx
        flags = [flags;1, 0];
    end
    if isy
        flags = [flags;0, 1];
    end
    if isx && isy
        flags = [flags;1, 1];
    end
    panelCount = size(vertices, 1);
    ni = size(flags, 1);
    all_vertices = zeros(panelCount, 4, 3, ni);
    all_normals = zeros(panelCount, 3, ni);
    weights = ones(ni, 1);
    labels = cell(ni, 1);
    for imageIndex = 1:ni
        imageVertices = vertices;
        normal = normals;
        if flags(imageIndex, 1)
            imageVertices(:, :, 1) = -imageVertices(:, :, 1);
            normal(:, 1) = -normal(:, 1);
        end
        if flags(imageIndex, 2)
            imageVertices(:, :, 2) = -imageVertices(:, :, 2);
            normal(:, 2) = -normal(:, 2);
        end
        if mod(sum(flags(imageIndex, :)), 2) == 1
            imageVertices(:, [2, 4], :) = imageVertices(:, [4, 2], :);
        end
        all_vertices(:, :, :, imageIndex) = imageVertices;
        all_normals(:, :, imageIndex) = normal;
        weights(imageIndex) = parity(1)^flags(imageIndex, 1) * parity(2)^flags(imageIndex, 2);
        labels{imageIndex} = sprintf('reflect_x=%d,reflect_y=%d', flags(imageIndex, 1), flags(imageIndex, 2));
    end
end
