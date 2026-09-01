function [rhs, sourceFlux, info] = build_rankine_radiation_rhs( ...
        singleLayerMatrix, omega, generalizedNormals, modeIndices)
% BUILD_RANKINE_RADIATION_RHS Build the radiation Green-identity RHS.
% The program convention is exp(+i*omega*t), and every radiation column
% represents unit generalized displacement.

    policy = radiation_rhs_policy();
    if nargin == 0
        rhs = complex(zeros(0, 0));
        sourceFlux = complex(zeros(0, 0));
        info = policy;
        return;
    end
    validateattributes(singleLayerMatrix, {'numeric'}, {'2d','finite'});
    validateattributes(omega, {'numeric'}, ...
        {'scalar','real','positive','finite'});
    validateattributes(generalizedNormals, {'numeric'}, {'2d','real','finite'});
    assert(size(singleLayerMatrix, 2) == size(generalizedNormals, 1), ...
        'CRESTU:RadiationRHSSize', ...
        'S source-column count must match generalized-normal rows.');
    if nargin < 4 || isempty(modeIndices)
        modeIndices = 1:size(generalizedNormals, 2);
    end
    modeIndices = reshape(modeIndices, 1, []);
    validateattributes(modeIndices, {'numeric'}, ...
        {'integer','positive','finite','numel',size(generalizedNormals, 2)});

    % n_s=n_B=-n_Omega.  With exp(+i*omega*t) and unit generalized
    % displacement, q_s=dphi/dn_s=i*omega*n_j.  Green's identity therefore
    % places the prescribed Neumann data on the right as -S*q_s.
    sourceFlux = 1i * omega * generalizedNormals;
    rhs = -singleLayerMatrix * sourceFlux;
    info = policy;
    info.omegaRadPerS = omega;
    info.equationCount = size(singleLayerMatrix, 1);
    info.sourcePanelCount = size(singleLayerMatrix, 2);
    info.modeIndices = modeIndices;
    info.modeCount = numel(modeIndices);
    info.sourceFluxSHA256 = numeric_sha256(sourceFlux);
    info.rightHandSideSHA256 = numeric_sha256(rhs);
end

function policy = radiation_rhs_policy()
    policy = struct( ...
        'modelName', 'radiation-green-identity-rhs-v1', ...
        'timeConvention', 'exp(+i*omega*t)', ...
        'motionNormalization', 'unit-generalized-displacement', ...
        'storedBodyNormalConvention', 'n_s=n_B=-n_Omega', ...
        'sourceFluxFormula', 'q_s=i*omega*n_j', ...
        'greenIdentityRHSFormula', 'rhs=-S*q_s', ...
        'diffractionRHSChanged', false, ...
        'postSolvePotentialSignAdjustment', false);
    policy.policySHA256 = sha256_hash(jsonencode(policy));
end

function hashText = numeric_sha256(values)
    values = double(values);
    signature = [double(size(values, 1)); double(size(values, 2)); ...
        real(values(:)); imag(values(:))];
    hashText = sha256_hash(typecast(signature, 'uint8'));
end
