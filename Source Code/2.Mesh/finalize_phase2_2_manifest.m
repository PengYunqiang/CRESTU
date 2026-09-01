function [cfg, manifest] = finalize_phase2_2_manifest(cfg, meshAudit, stats)
% FINALIZE_PHASE2_2_MANIFEST Bind actual component counts/hashes and radii.

    assert(isfield(cfg, 'phase2_2_controls') && ...
        isfield(cfg.phase2_2_controls, 'manifest'), ...
        'CRESTU:Phase22ManifestMissing', ...
        'Phase-2.2 controls must be resolved before manifest finalization.');
    manifest = cfg.phase2_2_controls.manifest;
    effective = manifest.effective;
    truncation = cfg.outer_truncation.effective;
    if isfield(truncation, 'effectiveTopRadiusM')
        topRadius = truncation.effectiveTopRadiusM;
        bottomRadius = truncation.effectiveBottomRadiusM;
    else
        topRadius = cfg.fs.r_outer;
        bottomRadius = cfg.fs.r_outer;
    end
    effective.topTruncationRadiusM = topRadius;
    effective.bottomTruncationRadiusM = bottomRadius;
    effective.spongePlateauWidthToTopM = max(0, ...
        topRadius - effective.spongePlateauStartRadiusM);
    effective.bodyPanelCount = stats.total_body_panels;
    effective.freeSurfacePanelCount = stats.fs_panels;
    effective.bottomPanelCount = stats.seabed_panels;
    effective.outerBoundaryPanelCount = stats.farfield_panels;
    effective.totalUnknownCount = stats.total_dofs;
    effective.bodyMeshSHA256 = meshAudit.body.hash;
    effective.freeSurfaceMeshSHA256 = meshAudit.freeSurface.hash;
    effective.bottomMeshSHA256 = meshAudit.bottom.hash;
    effective.outerBoundaryMeshSHA256 = meshAudit.outerBoundary.hash;
    effective.mergedGeometrySHA256 = meshAudit.geometry.hash;
    effective.outerTruncationInput = rmfield(cfg.outer_truncation, 'effective');
    effective.outerTruncationEffective = cfg.outer_truncation.effective;
    manifest.effective = effective;
    manifest.effectiveSHA256 = sha256_hash(jsonencode(effective));
    cfg.phase2_2_controls.effective = effective;
    cfg.phase2_2_controls.manifest = manifest;
end
