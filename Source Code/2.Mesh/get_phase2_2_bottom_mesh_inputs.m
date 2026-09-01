function [bottomCfg, bottomReferenceSurface] = ...
        get_phase2_2_bottom_mesh_inputs(cfg, waterline, freeSurface)
% GET_PHASE2_2_BOTTOM_MESH_INPUTS Decouple bottom radial counts for OFAT.
% The legacy-absent path returns the already generated free surface exactly.

    bottomCfg = cfg;
    bottomReferenceSurface = freeSurface;
    if ~isfield(cfg, 'phase2_2_controls') || ...
            ~cfg.phase2_2_controls.sectionPresent
        return
    end
    counts = cfg.phase2_2_controls.bottomRadialCounts;
    bottomCfg.fs.nr_near = counts(1);
    bottomCfg.fs.nr_sponge = counts(2);
    bottomReferenceSurface = generate_free_surface_bmf('', waterline, ...
        bottomCfg, false);
end
