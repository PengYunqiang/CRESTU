function expected = find_frozen_row(definition, level, omega)
% FIND_FROZEN_ROW Return exactly one v5 row for a level/frequency pair.

    frozenFile = fullfile(definition.artifactDirectory, ...
        'Phase3_2B_V5_Frozen_Mesh_Manifest.csv');
    frozen = readtable(frozenFile, 'TextType', 'string');
    frequencyTolerance = 64 * eps(max(1.0, abs(omega)));
    selected = frozen.mesh_level_index == level.index & ...
        abs(frozen.omega_rad_s - omega) <= frequencyTolerance;
    assert(nnz(selected) == 1, 'CRESTU:Phase32BFrozenRow', ...
        'Exactly one frozen row is required for %s at %.6g rad/s.', ...
        level.name, omega);
    expected = frozen(selected, :);
end
