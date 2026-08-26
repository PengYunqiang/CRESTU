function [velocity, hessian, diagnostics] = estimate_surface_kinematics(mesh, phi, dphi_dn, n_neighbors)
% ESTIMATE_SURFACE_KINEMATICS Recover surface velocity and potential Hessian by weighted local reconstruction.
%
% Syntax:
%   [velocity, hessian, diagnostics] = estimate_surface_kinematics(mesh, phi, dphi_dn, n_neighbors)
%
% Description:
%   The routine evaluates, reconstructs, imports, or exports quantities required by second-order mean wave-drift analysis. Complex products are time averaged consistently with the exp(i*omega*t) convention and generalized loads use the project 6-DOF ordering.
%
% Inputs:
%   mesh               - [struct] Boundary-panel mesh with Cartesian geometry in SI units.
%   phi                - [N x K] Complex velocity-potential samples, [m^2/s].
%   dphi_dn            - [N x K] Complex normal potential derivatives, [m/s].
%   n_neighbors        - [scalar] Number of panels in each local reconstruction stencil, dimensionless.
%
% Outputs:
%   velocity           - [N x 3 x K] Reconstructed complex surface velocity, [m/s].
%   hessian            - [N x 3 x 3 x K] Reconstructed Cartesian potential Hessian, [1/s].
%   diagnostics        - [struct] Numerical conditioning, symmetry, residual, or reconstruction diagnostics.
%
% Governing Equations / Theory:
%   Pinkster near-field direct pressure integration, Maruo-Newman far-field momentum balance, or supporting surface reconstruction as applicable.
%
% References:
%   - Pinkster, J. A. (1980), Low Frequency Second Order Wave Exciting Forces on Floating Structures; Newman, J. N. (1974), second-order slowly varying forces.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if nargin < 4 || isempty(n_neighbors)
        n_neighbors = 20;
    end
    panel_count = mesh.n_panels;
    heading_count = size(phi, 2);
    dphi_dn = reshape(dphi_dn, panel_count, heading_count);
    if size(phi, 1) ~= panel_count
        error('CRESTU:KinematicsShape', 'phi rows must equal the number of body panels.');
    end
    if panel_count < 6
        error('CRESTU:KinematicsStencil', 'At least six panels are required for surface reconstruction.');
    end

    centers = reshape(mesh.centers, panel_count, 3);
    normals = reshape(mesh.normals, panel_count, 3);
    tangent_one = reshape(mesh.e1, panel_count, 3);
    tangent_two = reshape(mesh.e2, panel_count, 3);
    neighbor_count = min(panel_count - 1, max(6, round(n_neighbors)));

%% --- 2. Build reusable local stencils ---

    distance_squared = sum(centers .^ 2, 2) + sum(centers .^ 2, 2).' - 2 * (centers * centers.');
    distance_squared = max(distance_squared, 0);
    [~, distance_order] = sort(distance_squared, 2, 'ascend');
    neighbors = distance_order(:, 2:(neighbor_count + 1));
    design_matrices = cell(panel_count, 1);
    stencil_weights = cell(panel_count, 1);
    condition_number = zeros(panel_count, 1);

    for panel_index = 1:panel_count
        neighbor_index = neighbors(panel_index, :);
        displacement = centers(neighbor_index, :) - centers(panel_index, :);
        local_u = displacement * tangent_one(panel_index, :).';
        local_v = displacement * tangent_two(panel_index, :).';
        radius = sqrt(local_u .^ 2 + local_v .^ 2);
        support_radius = max(radius(end), eps);
        weight = exp(-4 * (radius / support_radius) .^ 2);
        design = [local_u, local_v, 0.5 * local_u .^ 2, local_u .* local_v, 0.5 * local_v .^ 2];
        weighted_design = design .* sqrt(weight);
        design_matrices{panel_index} = design;
        stencil_weights{panel_index} = weight;
        singular_values = svd(weighted_design, 'econ');
        condition_number(panel_index) = singular_values(1) / max(singular_values(end), eps);
    end

%% --- 3. Recover gradients from a quadratic surface fit ---

    velocity = complex(zeros(panel_count, 3, heading_count));
    tangent_hessian = complex(zeros(panel_count, 3, 3, heading_count));
    fit_residual = zeros(panel_count, heading_count);
    regularization = 1.0e-10;

    for heading_index = 1:heading_count
        for panel_index = 1:panel_count
            neighbor_index = neighbors(panel_index, :);
            design = design_matrices{panel_index};
            weight = stencil_weights{panel_index};
            weighted_design = design .* sqrt(weight);
            potential_difference = phi(neighbor_index, heading_index) - phi(panel_index, heading_index);
            weighted_rhs = potential_difference .* sqrt(weight);
            normal_matrix = weighted_design' * weighted_design;
            scale = max(trace(real(normal_matrix)) / 5, 1);
            coefficient = (normal_matrix + regularization * scale * eye(5)) ...
                      \ (weighted_design' * weighted_rhs);

            velocity(panel_index, :, heading_index) = ...
                coefficient(1) * tangent_one(panel_index, :) ...
                + coefficient(2) * tangent_two(panel_index, :) ...
                + dphi_dn(panel_index, heading_index) * normals(panel_index, :);

            local_hessian = [coefficient(3), coefficient(4); coefficient(4), coefficient(5)];
            tangent_basis = [tangent_one(panel_index, :); tangent_two(panel_index, :)];
            tangent_hessian(panel_index, :, :, heading_index) = tangent_basis' * local_hessian * tangent_basis;
            fit_residual(panel_index, heading_index) = norm(weighted_design * coefficient - weighted_rhs) ...
                      / max(norm(weighted_rhs), eps);
        end
    end

%% --- 4. Recover the Cartesian Hessian from gradient variation ---

    hessian = complex(zeros(panel_count, 3, 3, heading_count));
    hessian_residual = zeros(panel_count, heading_count);
    for heading_index = 1:heading_count
        for panel_index = 1:panel_count
            neighbor_index = neighbors(panel_index, :);
            displacement = centers(neighbor_index, :) - centers(panel_index, :);
            radius = sqrt(sum(displacement .^ 2, 2));
            support_radius = max(radius(end), eps);
            weight = exp(-4 * (radius / support_radius) .^ 2);
            weighted_displacement = displacement .* sqrt(weight);
            velocity_difference = reshape(velocity(neighbor_index, :, heading_index), neighbor_count, 3) ...
                - velocity(panel_index, :, heading_index);
            weighted_difference = velocity_difference .* sqrt(weight);
            normal_matrix = weighted_displacement' * weighted_displacement;
            scale = max(trace(real(normal_matrix)) / 3, 1);
            gradient_matrix = (normal_matrix + regularization * scale * eye(3)) ...
                      \ (weighted_displacement' * weighted_difference);
            recovered_hessian = 0.5 * (gradient_matrix + gradient_matrix.');

            % Blend the direct tangent curvature with the ambient gradient-variation fit.
            direct_tangent = reshape(tangent_hessian(panel_index, :, :, heading_index), 3, 3);
            tangent_projector = eye(3) - normals(panel_index, :).' * normals(panel_index, :);
            recovered_hessian = recovered_hessian ...
                + tangent_projector * (direct_tangent - recovered_hessian) * tangent_projector;
            hessian(panel_index, :, :, heading_index) = 0.5 * (recovered_hessian + recovered_hessian.');
            hessian_residual(panel_index, heading_index) = ...
                norm(weighted_displacement * recovered_hessian.' - weighted_difference, 'fro') ...
                      / max(norm(weighted_difference, 'fro'), eps);
        end
    end

    diagnostics = struct( ...
        'neighbor_count', neighbor_count, ...
        'condition_number', condition_number, ...
        'median_condition_number', median(condition_number), ...
        'max_condition_number', max(condition_number), ...
        'potential_fit_residual', fit_residual, ...
        'hessian_fit_residual', hessian_residual);
end
