function report = run_validation_smoke()
%RUN_VALIDATION_SMOKE Validate portable setup and canonical geometry inputs.

    projectRoot = setup_crestu();
    configFile = fullfile(projectRoot, 'examples', 'single_sphere', ...
        'Case1_Smoke.cfg');
    domain = build_bmf_domain(configFile);
    assert(domain.stats.total_body_panels > 0, ...
        'CRESTU:SmokeEmptyBody', 'Smoke body mesh is empty.');
    assert(domain.stats.total_dofs > domain.stats.total_body_panels, ...
        'CRESTU:SmokeIncompleteDomain', ...
        'Smoke domain did not build free-surface/bottom/outer panels.');

    report = struct();
    report.projectRoot = projectRoot;
    report.configFile = configFile;
    report.bodyPanels = domain.stats.total_body_panels;
    report.freeSurfacePanels = domain.stats.fs_panels;
    report.seabedPanels = domain.stats.seabed_panels;
    report.farfieldPanels = domain.stats.farfield_panels;
    report.totalPanels = domain.stats.total_dofs;
    fprintf('[OK] Validation smoke passed | body=%d | total=%d panels\n', ...
        report.bodyPanels, report.totalPanels);
end
