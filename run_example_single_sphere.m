function results = run_example_single_sphere(mode)
%RUN_EXAMPLE_SINGLE_SPHERE Run the canonical release smoke example.
%   RESULTS = RUN_EXAMPLE_SINGLE_SPHERE() runs one frequency.
%   RESULTS = RUN_EXAMPLE_SINGLE_SPHERE('FULL') runs the 16-frequency Fine case.

    if nargin < 1 || isempty(mode)
        mode = 'SMOKE';
    end
    mode = upper(string(mode));
    assert(any(mode == ["SMOKE", "FULL"]), ...
        'CRESTU:ExampleMode', 'Mode must be SMOKE or FULL.');

    projectRoot = setup_crestu();
    exampleRoot = fullfile(projectRoot, 'examples', 'single_sphere');
    if mode == "SMOKE"
        configFile = fullfile(exampleRoot, 'Case1_Smoke.cfg');
    else
        configFile = fullfile(exampleRoot, 'Case1_Fine.cfg');
    end

    assert(isfile(configFile), 'CRESTU:ExampleConfigMissing', ...
        'Example configuration was not found: %s', configFile);
    fprintf('[INFO] Running CRESTU single-sphere example | mode=%s\n', mode);
    results = run_frequency_domain_case(configFile, ...
        skipPhysicalDiagnostics=true);
    fprintf('[OK] Example completed | frequencies=%d | panels=%d\n', ...
        numel(results.omegas), results.stats.total_dofs);
end
