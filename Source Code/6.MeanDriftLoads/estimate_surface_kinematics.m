function [velocity, hessian, diagnostics] = estimate_surface_kinematics(mesh, phi, dphi_dn, n_neighbors)
% ESTIMATE_SURFACE_KINEMATICS Recover the body-surface velocity and Hessian.
%
% Syntax:
%   [velocity, hessian, diagnostics] = estimate_surface_kinematics(mesh, phi, dphi_dn, n_neighbors)
%
% Inputs:
%   mesh        : [struct] Body mesh with panel centers, normals, and tangent bases, in m.
%   phi         : [N x Nh] Complex first-order velocity potential, in m^2/s.
%   dphi_dn     : [N x Nh] Prescribed normal derivative of the potential, in m/s.
%   n_neighbors : [scalar] Number of neighboring panels in each weighted least-squares stencil.
%
% Outputs:
%   velocity    : [N x 3 x Nh] Complex surface velocity, in m/s.
%   hessian     : [N x 3 x 3 x Nh] Symmetric Cartesian potential Hessian, in 1/s.
%   diagnostics : [struct] Stencil conditioning and reconstruction residuals.
%
% Mathematical Reference:
%   A weighted quadratic moving-least-squares surface reconstruction is used for tangential
%   derivatives. The Neumann boundary condition closes the normal component of grad(phi).

% ==========================================
% Validate and normalize the input arrays
% ==========================================
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

% ==========================================
% Build reusable local stencils
% ==========================================
    distance_squared = sum(centers.^2, 2) + sum(centers.^2, 2).' - 2 * (centers * centers.');
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
        radius = sqrt(local_u.^2 + local_v.^2);
        support_radius = max(radius(end), eps);
        weight = exp(-4 * (radius / support_radius).^2);
        design = [local_u, local_v, 0.5 * local_u.^2, local_u .* local_v, 0.5 * local_v.^2];
        weighted_design = design .* sqrt(weight);
        design_matrices{panel_index} = design;
        stencil_weights{panel_index} = weight;
        singular_values = svd(weighted_design, 'econ');
        condition_number(panel_index) = singular_values(1) / max(singular_values(end), eps);
    end

% ==========================================
% Recover gradients from a quadratic surface fit
% ==========================================
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

% ==========================================
% Recover the Cartesian Hessian from gradient variation
% ==========================================
    hessian = complex(zeros(panel_count, 3, 3, heading_count));
    hessian_residual = zeros(panel_count, heading_count);
    for heading_index = 1:heading_count
        for panel_index = 1:panel_count
            neighbor_index = neighbors(panel_index, :);
            displacement = centers(neighbor_index, :) - centers(panel_index, :);
            radius = sqrt(sum(displacement.^2, 2));
            support_radius = max(radius(end), eps);
            weight = exp(-4 * (radius / support_radius).^2);
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
