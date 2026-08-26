function drift = compute_drift_nearfield(mesh, waterline, state, cfg)
% COMPUTE_DRIFT_NEARFIELD Evaluate Pinkster-type mean drift loads by direct near-field integration.
%
% Syntax:
%   drift = compute_drift_nearfield(mesh, waterline, state, cfg)
%
% Description:
%   The routine evaluates, reconstructs, imports, or exports quantities required by second-order mean wave-drift analysis. Complex products are time averaged consistently with the exp(i*omega*t) convention and generalized loads use the project 6-DOF ordering.
%
% Inputs:
%   mesh               - [struct] Boundary-panel mesh with Cartesian geometry in SI units.
%   waterline          - [struct] Ordered waterline nodes and segment metadata, with coordinates in [m].
%   state              - [struct] Frequency-domain potential, motion, and load state in SI units.
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%
% Outputs:
%   drift              - [struct] Mean second-order generalized loads and audited component terms, [N] and [N m].
%
% Governing Equations / Theory:
%   Pinkster near-field direct pressure integration, Maruo-Newman far-field momentum balance, or supporting surface reconstruction as applicable.
%
% References:
%   - Pinkster, J. A. (1980), Low Frequency Second Order Wave Exciting Forces on Floating Structures; Newman, J. N. (1974), second-order slowly varying forces.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    required_fields = {'phi', 'dphi_dn', 'rao', 'first_order_force', 'omega', 'headings'};
    for field_index = 1:numel(required_fields)
        field_name = required_fields{field_index};
        if ~isfield(state, field_name)
            error('CRESTU:DriftState', 'Near-field state lacks %s.', field_name);
        end
    end

    panel_count = mesh.n_panels;
    heading_count = size(state.phi, 2);
    areas = reshape(mesh.areas, panel_count, 1);
    centers = reshape(mesh.centers, panel_count, 3);
    normals = reshape(mesh.normals, panel_count, 3);
    center_of_gravity = reshape(mesh.cg, 1, 3);
    generalized_normal = [normals, cross(centers - center_of_gravity, normals, 2)];

%% --- 2. Reconstruct smooth surface kinematics ---

    if ~isfield(state, 'velocity') || isempty(state.velocity) ...
                  || ~isfield(state, 'hessian') || isempty(state.hessian)
        [velocity, hessian, kinematics_diagnostics] = ...
            estimate_surface_kinematics(mesh, state.phi, state.dphi_dn, 20);
    else
        velocity = state.velocity;
        hessian = state.hessian;
        kinematics_diagnostics = struct('source', 'supplied');
    end

%% --- 3. Construct the exact outward waterline contour normal ---

    complete_waterline = waterline;
    if cfg.isx || cfg.isy
        complete_waterline = complete_waterline_by_symmetry(waterline, cfg.isx, cfg.isy);
    end
    nodes = complete_waterline.nodes;
    node_count = size(nodes, 1);
    next_node = [2:node_count, 1];
    segment = nodes(next_node, :) - nodes;
    segment_length = sqrt(sum(segment .^ 2, 2));
    if signed_polygon_area(nodes) < 0
        nodes = flipud(nodes);
        next_node = [2:node_count, 1];
        segment = nodes(next_node, :) - nodes;
        segment_length = sqrt(sum(segment .^ 2, 2));
    end
    midpoint = 0.5 * (nodes + nodes(next_node, :));
    line_normal = [segment(:, 2), -segment(:, 1)] ./ max(segment_length, eps);
    line_position = [midpoint, zeros(node_count, 1)] - center_of_gravity;
    line_normal_three = [line_normal, zeros(node_count, 1)];
    line_generalized_normal = [line_normal_three, cross(line_position, line_normal_three, 2)];

%% --- 4. Interpolate the potential continuously to waterline nodes ---

    if isfield(state, 'waterline_potential') && ~isempty(state.waterline_potential)
        waterline_potential = state.waterline_potential;
    else
        waterline_potential = interpolate_panel_values( ...
            [nodes, zeros(node_count, 1)], centers, state.phi, 6);
    end

%% --- 5. Evaluate the four time-averaged contributions ---

    term_waterline = zeros(6, heading_count);
    term_quadratic_velocity = zeros(6, heading_count);
    term_rotation_force = zeros(6, heading_count);
    term_translation_gradient = zeros(6, heading_count);

    for heading_index = 1:heading_count
        translation = state.rao(1:3, heading_index);
        rotation = state.rao(4:6, heading_index);

        wave_elevation = -1i * state.omega * waterline_potential(:, heading_index) / cfg.grav;
        body_elevation = translation(3) ...
            + rotation(1) * (nodes(:, 2) - center_of_gravity(2)) ...
            - rotation(2) * (nodes(:, 1) - center_of_gravity(1));
        relative_elevation = wave_elevation - body_elevation;
        midpoint_elevation = 0.5 * (relative_elevation + relative_elevation(next_node));
        term_waterline(:, heading_index) = 0.5 * cfg.rho * cfg.grav ...
                  * (line_generalized_normal.' * (abs(midpoint_elevation) .^ 2 .* segment_length));

        surface_velocity = reshape(velocity(:, :, heading_index), panel_count, 3);
        speed_squared = real(sum(surface_velocity .* conj(surface_velocity), 2));
        term_quadratic_velocity(:, heading_index) = -0.5 * cfg.rho ...
                  * (generalized_normal.' * (speed_squared .* areas));

        % The requested cross term uses the first-order excitation force only.
        term_rotation_force(1:3, heading_index) = 0.5 * real(cross( ...
            conj(rotation(:)), state.first_order_force(1:3, heading_index)));

        % Rotation is excluded here because its separate force cross-product is Term 3.
        translated_normal_acceleration = zeros(panel_count, 1);
        for panel_index = 1:panel_count
            potential_hessian = reshape(hessian(panel_index, :, :, heading_index), 3, 3);
            directional_gradient = potential_hessian * conj(translation(:));
            translated_normal_acceleration(panel_index) = real( ...
                dot(directional_gradient, normals(panel_index, :).'));
        end
        term_translation_gradient(:, heading_index) = -0.5 * cfg.rho ...
                  * (generalized_normal.' * (translated_normal_acceleration .* areas));
    end

    total = term_waterline + term_quadratic_velocity + term_rotation_force + term_translation_gradient;
    drift = struct( ...
        'method', 'near_field_direct', ...
        'omega', state.omega, ...
        'headings', state.headings, ...
        'term_waterline', term_waterline, ...
        'term_quadratic_velocity', term_quadratic_velocity, ...
        'term_rotation_force', term_rotation_force, ...
        'term_translation_gradient', term_translation_gradient, ...
        'total', total, ...
        'velocity', velocity, ...
        'hessian', hessian, ...
        'kinematics_diagnostics', kinematics_diagnostics, ...
        'waterline_orientation', 'counterclockwise with outward right-hand normal');
end

function area = signed_polygon_area(nodes)
% SIGNED_POLYGON_AREA Evaluate the signed area of an ordered planar polygon.
%
% Syntax:
%   area = signed_polygon_area(nodes)
%
% Description:
%   The routine evaluates, reconstructs, imports, or exports quantities required by second-order mean wave-drift analysis. Complex products are time averaged consistently with the exp(i*omega*t) convention and generalized loads use the project 6-DOF ordering.
%
% Inputs:
%   nodes              - [N x 2] Ordered planar polygon or waterline nodes, [m].
%
% Outputs:
%   area               - [scalar] Signed planar polygon area, [m^2].
%
% Governing Equations / Theory:
%   Pinkster near-field direct pressure integration, Maruo-Newman far-field momentum balance, or supporting surface reconstruction as applicable.
%
% References:
%   - Pinkster, J. A. (1980), Low Frequency Second Order Wave Exciting Forces on Floating Structures; Newman, J. N. (1974), second-order slowly varying forces.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    next_node = [2:size(nodes, 1), 1];
    area = 0.5 * sum(nodes(:, 1) .* nodes(next_node, 2) - nodes(next_node, 1) .* nodes(:, 2));
end

function values = interpolate_panel_values(points, centers, panel_values, neighbor_count)
% INTERPOLATE_PANEL_VALUES Interpolate panel-centered fields with inverse-distance weighting.
%
% Syntax:
%   values = interpolate_panel_values(points, centers, panel_values, neighbor_count)
%
% Description:
%   The routine evaluates, reconstructs, imports, or exports quantities required by second-order mean wave-drift analysis. Complex products are time averaged consistently with the exp(i*omega*t) convention and generalized loads use the project 6-DOF ordering.
%
% Inputs:
%   points             - [M x 3] Target interpolation points, [m].
%   centers            - [N x 3] Panel collocation points in global coordinates, [m].
%   panel_values       - [N x K] Panel-centered samples; physical units are preserved.
%   neighbor_count     - [scalar] Number of neighboring samples used for interpolation, dimensionless.
%
% Outputs:
%   values             - [M x K] Interpolated samples with the same physical units as the input field.
%
% Governing Equations / Theory:
%   Pinkster near-field direct pressure integration, Maruo-Newman far-field momentum balance, or supporting surface reconstruction as applicable.
%
% References:
%   - Pinkster, J. A. (1980), Low Frequency Second Order Wave Exciting Forces on Floating Structures; Newman, J. N. (1974), second-order slowly varying forces.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    distance_squared = sum(points .^ 2, 2) + sum(centers .^ 2, 2).' - 2 * (points * centers.');
    distance_squared = max(distance_squared, eps);
    [sorted_distance, sorted_index] = sort(distance_squared, 2, 'ascend');
    use_count = min(neighbor_count, size(centers, 1));
    sorted_distance = sorted_distance(:, 1:use_count);
    sorted_index = sorted_index(:, 1:use_count);
    weights = 1 ./ sorted_distance;
    weights = weights ./ sum(weights, 2);
    values = complex(zeros(size(points, 1), size(panel_values, 2)));
    for column_index = 1:size(panel_values, 2)
        local_values = reshape(panel_values(sorted_index, column_index), size(sorted_index));
        values(:, column_index) = sum(weights .* local_values, 2);
    end
end
