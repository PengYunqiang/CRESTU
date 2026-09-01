function report = audit_phase2_2_one_factor( ...
        baseline, candidate, expectedField, expectedOperatorEffect)
% AUDIT_PHASE2_2_ONE_FACTOR Verify one configured input and its whitelist.
% This helper reports failure rather than throwing so negative controls can
% prove that hidden co-variation is detected.

    if nargin < 4 || isempty(expectedOperatorEffect)
        expectedOperatorEffect = 'not-checked';
    end

    assert(isstruct(baseline) && isstruct(candidate) && ...
        isfield(baseline, 'independentInputs') && ...
        isfield(candidate, 'independentInputs'), ...
        'CRESTU:Phase22OFATManifest', ...
        'OFAT audit requires two configured-to-effective manifests.');
    expectedField = char(expectedField);
    baselineInputs = baseline.independentInputs;
    candidateInputs = candidate.independentInputs;
    names = fieldnames(baselineInputs);
    assert(isequal(names, fieldnames(candidateInputs)), ...
        'CRESTU:Phase22OFATShape', ...
        'Baseline and candidate independent-input fields differ.');
    changed = cell(0, 1);
    for nameIndex = 1:numel(names)
        name = names{nameIndex};
        if ~isequaln(baselineInputs.(name), candidateInputs.(name))
            changed{end + 1, 1} = name; %#ok<AGROW>
        end
    end
    inputPass = numel(changed) == 1 && strcmp(changed{1}, expectedField);
    [geometryPass, geometryChecks, rule] = check_geometry_whitelist( ...
        baseline, candidate, expectedField);
    [operatorPass, operatorChanged, operatorStatus, operatorHashes] = ...
        check_operator_effect(baseline, candidate, expectedOperatorEffect);
    report = struct('schemaVersion', 2, ...
        'expectedIndependentField', expectedField, ...
        'changedIndependentFields', {changed}, ...
        'independentInputPass', inputPass, ...
        'geometryWhitelistRule', rule, ...
        'geometryChecks', geometryChecks, ...
        'geometryWhitelistPass', geometryPass, ...
        'operatorExpectation', char(expectedOperatorEffect), ...
        'operatorChanged', operatorChanged, ...
        'operatorEffectStatus', operatorStatus, ...
        'operatorHashes', operatorHashes, ...
        'operatorPass', operatorPass, ...
        'passed', inputPass && geometryPass && operatorPass);
end

function [passed, changed, status, hashes] = check_operator_effect( ...
        baseline, candidate, expectation)
% CHECK_OPERATOR_EFFECT Distinguish configured levels from actual operators.
    field = 'boundaryOperatorStateSHA256';
    available = isfield(baseline, 'effective') && ...
        isfield(candidate, 'effective') && ...
        isfield(baseline.effective, field) && ...
        isfield(candidate.effective, field);
    hashes = struct('baseline', '', 'candidate', '');
    if ~available
        changed = false;
        status = 'NOT_AVAILABLE_UNBOUND_OPERATOR';
        passed = strcmpi(char(expectation), 'not-checked');
        return
    end
    hashes.baseline = baseline.effective.(field);
    hashes.candidate = candidate.effective.(field);
    changed = ~strcmp(hashes.baseline, hashes.candidate);
    switch lower(char(expectation))
        case 'changed'
            passed = changed;
            if changed
                status = 'ACTUAL_OPERATOR_CHANGED';
            else
                status = 'NO_EFFECT';
            end
        case 'unchanged'
            passed = ~changed;
            if changed
                status = 'UNEXPECTED_OPERATOR_CHANGE';
            else
                status = 'NO_EFFECT_EXPECTED';
            end
        case 'scientific-level'
            passed = changed;
            if changed
                status = 'ACTUAL_OPERATOR_CHANGED';
            else
                status = 'NO_EFFECT/NOT_A_SCIENTIFIC_LEVEL';
            end
        case 'not-checked'
            passed = true;
            status = 'NOT_CHECKED';
        otherwise
            error('CRESTU:Phase22OperatorExpectation', ...
                'Unknown operator-effect expectation %s.', ...
                char(expectation));
    end
end

function [passed, checks, rule] = check_geometry_whitelist( ...
        baseline, candidate, expectedField)
    required = {'bodyMeshSHA256', 'freeSurfaceMeshSHA256', ...
        'bottomMeshSHA256', 'outerBoundaryMeshSHA256', ...
        'mergedGeometrySHA256'};
    if ~isfield(baseline, 'effective') || ~isfield(candidate, 'effective') || ...
            ~all(cellfun(@(name) isfield(baseline.effective, name) && ...
            isfield(candidate.effective, name), required))
        checks = struct('status', 'NOT_AVAILABLE_UNFINALIZED');
        rule = 'finalized component hashes required';
        passed = false;
        return
    end
    first = baseline.effective;
    second = candidate.effective;
    sameBody = strcmp(first.bodyMeshSHA256, second.bodyMeshSHA256);
    sameFS = strcmp(first.freeSurfaceMeshSHA256, second.freeSurfaceMeshSHA256);
    sameBottom = strcmp(first.bottomMeshSHA256, second.bottomMeshSHA256);
    sameOuter = strcmp(first.outerBoundaryMeshSHA256, second.outerBoundaryMeshSHA256);
    sameMerged = strcmp(first.mergedGeometrySHA256, second.mergedGeometrySHA256);
    checks = struct('sameBodyHash', sameBody, ...
        'sameFreeSurfaceHash', sameFS, 'sameBottomHash', sameBottom, ...
        'sameOuterBoundaryHash', sameOuter, ...
        'sameMergedGeometryHash', sameMerged);
    switch expectedField
        case 'fsRadialCounts'
            rule = '#1 FS spacing: only FS component geometry may change';
            passed = sameBody && ~sameFS && sameBottom && sameOuter && ~sameMerged;
        case 'bottomRadialCounts'
            rule = '#3 bottom spacing: only bottom component geometry may change';
            passed = sameBody && sameFS && ~sameBottom && sameOuter && ~sameMerged;
        case 'outerVerticalPanelCount'
            rule = '#5 outer vertical spacing: only outer component geometry may change';
            passed = sameBody && sameFS && sameBottom && ~sameOuter && ~sameMerged;
        case 'outerThetaPanelCount'
            rule = 'common endpoint theta resolution: FS, bottom, and outer may change';
            passed = sameBody && ~sameFS && ~sameBottom && ~sameOuter && ~sameMerged;
        case 'freeSurfaceOuterThetaPanelCount'
            rule = 'FS endpoint theta resolution: FS and outer may change; bottom fixed';
            passed = sameBody && ~sameFS && sameBottom && ~sameOuter && ~sameMerged;
        case 'bottomOuterThetaPanelCount'
            rule = 'bottom endpoint theta resolution: bottom and outer may change; FS fixed';
            passed = sameBody && sameFS && ~sameBottom && ~sameOuter && ~sameMerged;
        case 'meshTransitionRadiusM'
            rule = 'transition radius may change FS and bottom, never body/outer';
            passed = sameBody && ~sameFS && ~sameBottom && sameOuter && ~sameMerged;
        otherwise
            rule = 'BC/diagnostic factor: all component and merged geometry hashes fixed';
            passed = sameBody && sameFS && sameBottom && sameOuter && sameMerged;
    end
end
