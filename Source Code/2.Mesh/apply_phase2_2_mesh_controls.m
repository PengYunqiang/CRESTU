function cfg = apply_phase2_2_mesh_controls(cfg)
% APPLY_PHASE2_2_MESH_CONTROLS Apply optional, geometry-only audit controls.
% The legacy-absent path returns cfg without changing any legacy fs value.

    if ~isfield(cfg, 'phase2_2_controls') || ...
            ~cfg.phase2_2_controls.sectionPresent
        return
    end
    controls = cfg.phase2_2_controls;
    assert(cfg.n_bodies == 1, 'CRESTU:Phase22ControlScope', ...
        'Active Phase-2.2 domain controls are limited to the single-body audit.');
    assert(controls.meshTransitionRadiusM > 0 && ...
        controls.meshTransitionRadiusM < cfg.fs.r_outer, ...
        'CRESTU:Phase22TransitionRadius', ...
        'MESH_TRANSITION_RADIUS must lie strictly inside the active outer radius.');
    cfg.fs.nr_near = controls.fsRadialCounts(1);
    cfg.fs.nr_sponge = controls.fsRadialCounts(2);
    cfg.fs.nz_farfield = controls.outerVerticalPanelCount;
end
