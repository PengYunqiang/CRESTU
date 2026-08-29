function rao = solve_rao(omegas, added_mass, damping, hydrostatic, excitation, cfg)
% SOLVE_RAO Solve the coupled frequency-domain rigid-body equations for response amplitude operators.
%
% Syntax:
%   rao = solve_rao(omegas, added_mass, damping, hydrostatic, excitation, cfg)
%
% Description:
%   Computes hydrodynamic coefficients, loads, or rigid-body response.
%   Results retain the CRESTU global 6-DOF order.
%
% Inputs:
%   omegas             - [1 x Nf] Angular-frequency samples, [rad/s].
%   added_mass         - [Ndof x Ndof x Nf] Added-mass coefficients in consistent translational and rotational SI units.
%   damping            - [Ndof x Ndof x Nf] Radiation-damping coefficients in consistent translational and rotational SI units.
%   hydrostatic        - [Ndof x Ndof] Linear hydrostatic restoring matrix in consistent SI units.
%   excitation         - [Ndof x Nh x Nf] Complex first-order wave-excitation loads, [N] and [N m].
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%
% Outputs:
%   rao                - [Ndof x Nh x Nf] Complex motion response per unit wave amplitude, [m/m] and [rad/m].
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    omegas = reshape(omegas, 1, []);
    frequencyCount = numel(omegas);
    globalDofCount = 6 * cfg.n_bodies;
    assert(frequencyCount > 0 && all(omegas > 0), ...
        'CRESTU:RAOFrequency', 'All angular frequencies must be positive.');
    added_mass = normalize_coeff(added_mass, globalDofCount, frequencyCount, 'added mass');
    damping = normalize_coeff(damping, globalDofCount, frequencyCount, 'damping');

    % A five-module model has 30 global DOFs. Module m uses
    % 6*(m-1)+(1:6). Local DOF 1-6: Surge, Sway, Heave, Roll, Pitch, Yaw.
    if ~isequal(size(hydrostatic), [globalDofCount, globalDofCount])
        error('CRESTU:HydrostaticShape', ...
            'Hydrostatic matrix must be %d-by-%d.', globalDofCount, globalDofCount);
    end

    if ismatrix(excitation) && frequencyCount == 1
        excitation = reshape(excitation, globalDofCount, size(excitation, 2), 1);
    end

    if size(excitation, 1) ~= globalDofCount || ...
            size(excitation, 3) ~= frequencyCount
        error('CRESTU:ExcitationShape', ...
            'Excitation must be ndof-by-nheading-by-nfrequency.');
    end

    %% Stage 2: Solve the coupled dynamic equation at each frequency

    headingCount = size(excitation, 2);
    M = assemble_mass_matrix(cfg.mass_props);
    complexResponse = complex(zeros(globalDofCount, headingCount, frequencyCount));
    dynamicMatrix = complex(zeros(globalDofCount, globalDofCount, frequencyCount));
    reciprocalCondition = zeros(frequencyCount, 1);

    % <<<CORE>>> solve_frequency_domain_rao, paper_eq=rigid_body_frequency_equation, benchmark=single_sphere_rao
    for k = 1:frequencyCount
        omega = omegas(k); % [rad/s]
        dynamicStiffness = -omega ^ 2 * (M + added_mass(:, :, k)) + ...
            1i * omega * damping(:, :, k) + hydrostatic;
        dynamicMatrix(:, :, k) = dynamicStiffness;
        reciprocalCondition(k) = rcond(dynamicStiffness);

        if reciprocalCondition(k) < 1e-12
            warning('CRESTU:IllConditionedRAO', ...
                'RAO matrix at omega=%g has rcond=%g.', ...
                omega, reciprocalCondition(k));
        end

        complexResponse(:, :, k) = dynamicStiffness \ excitation(:, :, k);
    end
    % <<</CORE>>>

    %% Stage 3: Preserve the response-field contract

    rao = struct('complex', complexResponse, 'amplitude', abs(complexResponse), ...
        'phase_deg', angle(complexResponse) * 180 / pi, ...
        'mass_matrix', M, 'dynamic_matrix', dynamicMatrix, ...
        'rcond', reciprocalCondition, 'omegas', omegas, ...
        'headings', cfg.wave.headings);
end

function out = normalize_coeff(in, n, nf, label)
% NORMALIZE_COEFF Normalize frequency-dependent coefficient arrays to the required shape.
%
% Syntax:
%   out = normalize_coeff(in, n, nf, label)
%
% Description:
%   Computes hydrodynamic coefficients, loads, or rigid-body response.
%   Results retain the CRESTU global 6-DOF order.
%
% Inputs:
%   in                 - [character vector, string scalar, or numeric value] Value to normalize.
%   n                  - [scalar] Number of coupled degrees of freedom or array rows, dimensionless.
%   nf                 - [scalar] Number of frequency samples, dimensionless.
%   label              - [character vector or string scalar] Diagnostic field label.
%
% Outputs:
%   out                - [value] Normalized or selected result with type and physical units preserved from the input.
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if iscell(in)
        out = zeros(n, n, nf);
        for k = 1:nf
            out(:, :, k) = in{k};
        end
    else
        out = in;
    end
    if nf == 1 && ismatrix(out)
        out = reshape(out, n, n, 1);
    end
    if size(out, 1) ~= n || size(out, 2) ~= n || size(out, 3) ~= nf
        error('CRESTU:CoefficientShape','%s array has the wrong shape.', label);
    end
end
