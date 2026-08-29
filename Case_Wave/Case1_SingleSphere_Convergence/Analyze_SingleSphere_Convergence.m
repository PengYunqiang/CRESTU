function study = Analyze_SingleSphere_Convergence( ...
        studyName, studyType, levelNames, resultSet, outputDirectory)
% ANALYZE_SINGLESPHERE_CONVERGENCE Quantify real body or outer-mesh convergence.
%
% Syntax:
%   study = Analyze_SingleSphere_Convergence(...)
%
% Inputs:
%   studyName      - Filename/title prefix [-].
%   studyType      - 'body' or 'outer' characteristic-size selector [-].
%   levelNames     - Coarse/Medium/Fine labels [-].
%   resultSet      - Three first-order CRESTU result structures [-].
%   outputDirectory - Directory for CSV, MAT, and PNG outputs [-].
%
% Outputs:
%   study          - Curves, differences, extrapolates, errors, and checks.

    arguments
        studyName {mustBeTextScalar}
        studyType {mustBeTextScalar}
        levelNames cell
        resultSet cell
        outputDirectory {mustBeTextScalar}
    end

    %% 阶段 1: 验证三层结果并提取 A33/B33/F3/RAO

    assert(numel(levelNames) == 3 && numel(resultSet) == 3, ...
        'CRESTU:ConvergenceLevelCount', ...
        'A convergence study requires exactly Coarse/Medium/Fine results.');
    assert(isfolder(outputDirectory), 'CRESTU:ConvergenceOutputDirectory', ...
        'Convergence output directory does not exist: %s', outputDirectory);
    omega = reshape(resultSet{1}.omegas, [], 1); % [rad/s]
    frequencyCount = numel(omega);
    quantityNames = {'A33', 'B33', 'F3', 'RAO3'};
    quantityUnits = {'kg', 'kg/s', 'N', 'm/m'};
    curves = zeros(frequencyCount, 3, 4);
    characteristicSize = zeros(frequencyCount, 3); % [m]

    for levelIndex = 1:3
        result = resultSet{levelIndex};
        assert(max(abs(reshape(result.omegas, [], 1) - omega)) < 1.0e-12, ...
            'CRESTU:ConvergenceFrequencyGrid', ...
            'All convergence levels must use the same frequency grid.');
        headingIndex = find(abs(result.headings) < 1.0e-12, 1);
        assert(~isempty(headingIndex), 'CRESTU:ConvergenceHeading', ...
            'Zero-degree wave heading is required for F3 and RAO3.');
        curves(:, levelIndex, 1) = squeeze(result.added_mass(3, 3, :));
        curves(:, levelIndex, 2) = squeeze(result.damping(3, 3, :));
        curves(:, levelIndex, 3) = squeeze(abs( ...
            result.excitation(3, headingIndex, :)));
        curves(:, levelIndex, 4) = squeeze( ...
            result.rao.amplitude(3, headingIndex, :));
        characteristicSize(:, levelIndex) = select_mesh_size( ...
            result.audit.frequencyEntries, studyType);
    end

    %% 阶段 2: 计算层间差、零网格尺度外推与自动判据

    coarseMediumDifference = squeeze(curves(:, 1, :) - curves(:, 2, :));
    mediumFineDifference = squeeze(curves(:, 2, :) - curves(:, 3, :));
    extrapolated = zeros(frequencyCount, 4);
    absoluteError = zeros(frequencyCount, 3, 4);

    for frequencyIndex = 1:frequencyCount
        meshSizes = characteristicSize(frequencyIndex, :); % [m]

        for quantityIndex = 1:4
            quantityValues = curves(frequencyIndex, :, quantityIndex);
            extrapolated(frequencyIndex, quantityIndex) = ...
                extrapolate_zero_mesh_size(meshSizes, quantityValues);
            absoluteError(frequencyIndex, :, quantityIndex) = abs( ...
                quantityValues - extrapolated(frequencyIndex, quantityIndex));
        end
    end

    checkTable = build_check_table(quantityNames, curves, absoluteError, ...
        coarseMediumDifference, mediumFineDifference);
    longTable = build_long_table(omega, levelNames, quantityNames, ...
        characteristicSize, curves, extrapolated, absoluteError, ...
        coarseMediumDifference, mediumFineDifference);

    %% 阶段 3: 导出曲线、差值、网格尺度图和结果文件

    prefix = char(studyName);
    longFilename = fullfile(outputDirectory, ...
        sprintf('%s_Convergence_Long.csv', prefix));
    checkFilename = fullfile(outputDirectory, ...
        sprintf('%s_Convergence_Checks.csv', prefix));
    writetable(longTable, longFilename);
    writetable(checkTable, checkFilename);
    figureFiles = create_convergence_figures(prefix, omega, levelNames, ...
        quantityNames, quantityUnits, characteristicSize, curves, ...
        extrapolated, absoluteError, coarseMediumDifference, ...
        mediumFineDifference, outputDirectory);
    study = struct('name', prefix, 'type', char(studyType), ...
        'levels', {levelNames}, 'omegas', omega, ...
        'quantityNames', {quantityNames}, 'quantityUnits', {quantityUnits}, ...
        'characteristicSizeM', characteristicSize, 'curves', curves, ...
        'coarseMediumDifference', coarseMediumDifference, ...
        'mediumFineDifference', mediumFineDifference, ...
        'extrapolatedZeroMeshSize', extrapolated, ...
        'absoluteErrorToExtrapolate', absoluteError, ...
        'checks', checkTable, 'longTableFile', longFilename, ...
        'checkTableFile', checkFilename, 'figureFiles', {figureFiles});
    studyFilename = fullfile(outputDirectory, sprintf('%s_Study.mat', prefix));
    save(studyFilename, 'study', '-v7.3');
    fprintf('[INFO] %s convergence checks:\n', prefix);
    disp(checkTable);
    fprintf('[OK] %s convergence study saved: %s\n', prefix, studyFilename);
end

function meshSize = select_mesh_size(auditEntries, studyType)
% SELECT_MESH_SIZE Extract the characteristic size used by this study.

    if strcmpi(studyType, 'body')
        meshSize = reshape([auditEntries.bodyCharacteristicSizeM], [], 1); % [m]
    elseif strcmpi(studyType, 'outer')
        meshSize = reshape([auditEntries.outerCharacteristicSizeM], [], 1); % [m]
    else
        error('CRESTU:ConvergenceStudyType', ...
            'Study type must be body or outer.');
    end
end

function extrapolatedValue = extrapolate_zero_mesh_size(meshSizes, quantityValues)
% EXTRAPOLATE_ZERO_MESH_SIZE Fit q(h)=q0+c1*h+c2*h^2 through three levels.

    assert(numel(unique(meshSizes)) == 3, ...
        'CRESTU:ConvergenceMeshSizes', ...
        'The three characteristic mesh sizes must be distinct.');
    polynomialCoefficients = polyfit(meshSizes, quantityValues, 2);
    extrapolatedValue = polyval(polynomialCoefficients, 0.0);
end

function checkTable = build_check_table(quantityNames, curves, ...
        absoluteError, coarseMediumDifference, mediumFineDifference)
% BUILD_CHECK_TABLE Test whether Fine is closer to the extrapolate than Medium.

    quantityCount = numel(quantityNames);
    medianMediumError = zeros(quantityCount, 1);
    medianFineError = zeros(quantityCount, 1);
    rmsMediumError = zeros(quantityCount, 1);
    rmsFineError = zeros(quantityCount, 1);
    medianDifferenceRatio = zeros(quantityCount, 1);
    fineCloserFrequencyFraction = zeros(quantityCount, 1);
    fineCloserOverall = false(quantityCount, 1);

    for quantityIndex = 1:quantityCount
        mediumError = absoluteError(:, 2, quantityIndex);
        fineError = absoluteError(:, 3, quantityIndex);
        mediumScale = max(abs(curves(:, 2, quantityIndex)), eps);
        fineScale = max(abs(curves(:, 3, quantityIndex)), eps);
        relativeMediumError = mediumError ./ mediumScale;
        relativeFineError = fineError ./ fineScale;
        medianMediumError(quantityIndex) = median(relativeMediumError);
        medianFineError(quantityIndex) = median(relativeFineError);
        rmsMediumError(quantityIndex) = sqrt(mean(relativeMediumError.^2));
        rmsFineError(quantityIndex) = sqrt(mean(relativeFineError.^2));
        medianDifferenceRatio(quantityIndex) = median( ...
            abs(mediumFineDifference(:, quantityIndex)) ./ ...
            max(abs(coarseMediumDifference(:, quantityIndex)), eps));
        fineCloserFrequencyFraction(quantityIndex) = mean( ...
            fineError <= mediumError);
        fineCloserOverall(quantityIndex) = ...
            rmsFineError(quantityIndex) <= rmsMediumError(quantityIndex);
    end

    checkTable = table(string(quantityNames(:)), medianMediumError, ...
        medianFineError, rmsMediumError, rmsFineError, ...
        medianDifferenceRatio, fineCloserFrequencyFraction, ...
        fineCloserOverall, 'VariableNames', {'quantity', ...
        'medianMediumRelativeError', 'medianFineRelativeError', ...
        'rmsMediumRelativeError', 'rmsFineRelativeError', ...
        'medianAbsMediumFineOverCoarseMedium', ...
        'fineCloserFrequencyFraction', 'fineCloserOverall'});
end

function longTable = build_long_table(omega, levelNames, quantityNames, ...
        characteristicSize, curves, extrapolated, absoluteError, ...
        coarseMediumDifference, mediumFineDifference)
% BUILD_LONG_TABLE Export every frequency, quantity, and grid level.

    frequencyCount = numel(omega);
    rowCount = frequencyCount * 3 * 4;
    studyOmega = zeros(rowCount, 1); % [rad/s]
    level = strings(rowCount, 1);
    quantity = strings(rowCount, 1);
    meshSize = zeros(rowCount, 1); % [m]
    value = zeros(rowCount, 1);
    zeroMeshExtrapolate = zeros(rowCount, 1);
    absoluteErrorToExtrapolate = zeros(rowCount, 1);
    coarseMinusMedium = zeros(rowCount, 1);
    mediumMinusFine = zeros(rowCount, 1);
    rowIndex = 0;

    for frequencyIndex = 1:frequencyCount
        for levelIndex = 1:3
            for quantityIndex = 1:4
                rowIndex = rowIndex + 1;
                studyOmega(rowIndex) = omega(frequencyIndex);
                level(rowIndex) = string(levelNames{levelIndex});
                quantity(rowIndex) = string(quantityNames{quantityIndex});
                meshSize(rowIndex) = characteristicSize(frequencyIndex, levelIndex);
                value(rowIndex) = curves(frequencyIndex, levelIndex, quantityIndex);
                zeroMeshExtrapolate(rowIndex) = ...
                    extrapolated(frequencyIndex, quantityIndex);
                absoluteErrorToExtrapolate(rowIndex) = ...
                    absoluteError(frequencyIndex, levelIndex, quantityIndex);
                coarseMinusMedium(rowIndex) = ...
                    coarseMediumDifference(frequencyIndex, quantityIndex);
                mediumMinusFine(rowIndex) = ...
                    mediumFineDifference(frequencyIndex, quantityIndex);
            end
        end
    end

    longTable = table(studyOmega, level, quantity, meshSize, value, ...
        zeroMeshExtrapolate, absoluteErrorToExtrapolate, ...
        coarseMinusMedium, mediumMinusFine, 'VariableNames', ...
        {'omegaRadPerSecond', 'level', 'quantity', ...
        'characteristicMeshSizeM', 'value', 'zeroMeshExtrapolate', ...
        'absoluteErrorToExtrapolate', 'coarseMinusMedium', ...
        'mediumMinusFine'});
end

function figureFiles = create_convergence_figures(prefix, omega, ...
        levelNames, quantityNames, quantityUnits, characteristicSize, ...
        curves, extrapolated, absoluteError, coarseMediumDifference, ...
        mediumFineDifference, outputDirectory)
% CREATE_CONVERGENCE_FIGURES Create frequency, difference, and h-space plots.

    frequencyFigure = figure('Color', 'w', 'Visible', 'off', ...
        'Position', [80, 80, 1400, 900]);
    frequencyLayout = tiledlayout(frequencyFigure, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    colors = lines(3);

    for quantityIndex = 1:4
        nexttile(frequencyLayout);
        hold on

        for levelIndex = 1:3
            plot(omega, curves(:, levelIndex, quantityIndex), '-o', ...
                'Color', colors(levelIndex, :), ...
                'DisplayName', levelNames{levelIndex});
        end

        format_axes(quantityNames{quantityIndex}, quantityUnits{quantityIndex});
    end

    title(frequencyLayout, sprintf('%s: three-grid frequency curves', prefix));
    frequencyFilename = fullfile(outputDirectory, ...
        sprintf('%s_Frequency_Curves.png', prefix));
    exportgraphics(frequencyFigure, frequencyFilename, 'Resolution', 180);
    close(frequencyFigure);

    differenceFigure = figure('Color', 'w', 'Visible', 'off', ...
        'Position', [80, 80, 1400, 900]);
    differenceLayout = tiledlayout(differenceFigure, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    for quantityIndex = 1:4
        nexttile(differenceLayout);
        plot(omega, coarseMediumDifference(:, quantityIndex), '-o', ...
            'DisplayName', 'Coarse - Medium');
        hold on
        plot(omega, mediumFineDifference(:, quantityIndex), '-s', ...
            'DisplayName', 'Medium - Fine');
        yline(0.0, 'k:');
        format_axes([quantityNames{quantityIndex}, ' differences'], ...
            quantityUnits{quantityIndex});
    end

    title(differenceLayout, sprintf('%s: adjacent-grid differences', prefix));
    differenceFilename = fullfile(outputDirectory, ...
        sprintf('%s_Grid_Differences.png', prefix));
    exportgraphics(differenceFigure, differenceFilename, 'Resolution', 180);
    close(differenceFigure);

    meshFigure = figure('Color', 'w', 'Visible', 'off', ...
        'Position', [80, 80, 1400, 900]);
    meshLayout = tiledlayout(meshFigure, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    finalFrequencyIndex = numel(omega);

    for quantityIndex = 1:4
        nexttile(meshLayout);
        meshSizes = characteristicSize(finalFrequencyIndex, :); % [m]
        [meshSizes, sortIndex] = sort(meshSizes);
        quantityValues = curves(finalFrequencyIndex, sortIndex, quantityIndex);
        relativeErrors = absoluteError(finalFrequencyIndex, sortIndex, ...
            quantityIndex) / max(abs(extrapolated(finalFrequencyIndex, ...
            quantityIndex)), eps);
        yyaxis left
        plot(meshSizes, quantityValues, '-o', 'LineWidth', 1.2);
        ylabel(sprintf('%s (%s)', quantityNames{quantityIndex}, ...
            quantityUnits{quantityIndex}));
        yyaxis right
        semilogy(meshSizes, max(relativeErrors, eps), '--s', ...
            'LineWidth', 1.2);
        ylabel('relative error to h=0 extrapolate');
        xlabel('characteristic mesh size h (m)');
        title(sprintf('%s at \\omega=%.2f rad/s', ...
            quantityNames{quantityIndex}, omega(finalFrequencyIndex)));
        grid on
        box on
    end

    title(meshLayout, sprintf('%s: quantity/error versus mesh size', prefix));
    meshFilename = fullfile(outputDirectory, ...
        sprintf('%s_h_Convergence.png', prefix));
    exportgraphics(meshFigure, meshFilename, 'Resolution', 180);
    close(meshFigure);
    figureFiles = {frequencyFilename, differenceFilename, meshFilename};
end

function format_axes(plotTitle, quantityUnit)
% FORMAT_AXES Apply consistent formatting to convergence axes.

    grid on
    box on
    xlabel('\omega (rad/s)');
    ylabel(quantityUnit);
    title(plotTitle);
    legend('Location', 'best');
end
