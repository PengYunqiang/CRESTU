function auditTable = audit_code_quality(projectRoot, outputFile)
% AUDIT_CODE_QUALITY Audit MATLAB readability rules and export a CSV report.
%
% Syntax:
%   auditTable = audit_code_quality(projectRoot, outputFile)
%
% Description:
%   Scans MATLAB files inside one authorized project root. The audit reports
%   stage sections, compact control statements, non-ASCII text, and long lines.
%
% Inputs:
%   projectRoot - Authorized project root, character vector or string scalar [-].
%   outputFile  - Reconciliation CSV path, character vector or string scalar [-].
%
% Outputs:
%   auditTable  - One-row-per-file code-quality reconciliation table [-].
%
% Governing Equations / Theory:
%   This utility applies deterministic text metrics and does not change code.
%
% References:
%   - CRESTU engineering code standard.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

    %% Stage 1: Apply defaults and validate directory boundaries

    if nargin < 1 || isempty(projectRoot)
        toolDirectory = fileparts(mfilename('fullpath'));
        sourceDirectory = fileparts(toolDirectory);
        projectRoot = fileparts(sourceDirectory);
    end

    if nargin < 2 || isempty(outputFile)
        outputFile = fullfile(projectRoot, 'MATLAB_Code_Quality_Audit.csv');
    end

    projectRoot = char(projectRoot);
    outputFile = char(outputFile);
    assert(isfolder(projectRoot), 'CRESTU:AuditProjectRoot', ...
        'Project root was not found: %s', projectRoot);
    resolvedRoot = char(java.io.File(projectRoot).getCanonicalPath());
    resolvedOutput = char(java.io.File(outputFile).getCanonicalPath());
    assert(startsWith(resolvedOutput, resolvedRoot, 'IgnoreCase', ispc), ...
        'CRESTU:AuditOutputBoundary', ...
        'Audit output must remain inside the project root.');

    %% Stage 2: Scan MATLAB source files

    matlabFiles = dir(fullfile(resolvedRoot, '**', '*.m'));
    fileCount = numel(matlabFiles);
    relativeFile = strings(fileCount, 1);
    lineCount = zeros(fileCount, 1);
    functionCount = zeros(fileCount, 1);
    stageCount = zeros(fileCount, 1);
    compactControlCount = zeros(fileCount, 1);
    unclearSingleLetterCount = zeros(fileCount, 1);
    nonAsciiLineCount = zeros(fileCount, 1);
    longLineCount = zeros(fileCount, 1);

    fprintf('[INFO] Audit MATLAB code quality | files = %d\n', fileCount);

    for fileIndex = 1:fileCount
        absoluteFile = fullfile(matlabFiles(fileIndex).folder, ...
            matlabFiles(fileIndex).name);
        sourceText = fileread(absoluteFile);
        sourceLines = splitlines(string(sourceText));
        relativeFile(fileIndex) = erase(string(absoluteFile), ...
            string(resolvedRoot) + filesep);
        lineCount(fileIndex) = numel(sourceLines);
        functionCount(fileIndex) = numel(regexp(sourceText, ...
            '(?m)^\s*function\b', 'match'));
        stageCount(fileIndex) = count(sourceText, '%% Stage ');
        compactControlCount(fileIndex) = count_compact_controls(sourceLines);
        unclearSingleLetterCount(fileIndex) = count_unclear_single_letters(sourceLines);
        nonAsciiLineCount(fileIndex) = count_non_ascii_lines(sourceLines);
        longLineCount(fileIndex) = sum(strlength(sourceLines) > 120);
    end

    %% Stage 3: Export and print the reconciliation table

    status = repmat("OK", fileCount, 1);
    status(compactControlCount > 0 | unclearSingleLetterCount > 0 | ...
        nonAsciiLineCount > 0 | stageCount == 0) = "CHECK";
    auditTable = table(relativeFile, lineCount, functionCount, stageCount, ...
        compactControlCount, unclearSingleLetterCount, nonAsciiLineCount, ...
        longLineCount, status, ...
        'VariableNames', {'File', 'Lines', 'Functions', 'Stages', ...
        'CompactControls', 'UnclearSingleLetters', 'NonAsciiLines', ...
        'LongLines', 'Status'});
    writetable(auditTable, resolvedOutput);

    summaryTable = table(fileCount, sum(lineCount), sum(functionCount), ...
        sum(stageCount), sum(compactControlCount), sum(unclearSingleLetterCount), ...
        sum(nonAsciiLineCount), sum(longLineCount), ...
        'VariableNames', {'Files', 'Lines', 'Functions', 'Stages', ...
        'CompactControls', 'UnclearSingleLetters', 'NonAsciiLines', 'LongLines'});
    disp(summaryTable);

    if any(status == "CHECK")
        fprintf('[WARN] Code-quality audit completed with items to review.\n');
    else
        fprintf('[OK] Code-quality audit completed with no blocking items.\n');
    end

    fprintf('[OK] Audit CSV exported | file = %s\n', resolvedOutput);
end

function unclearSingleLetterCount = count_unclear_single_letters(sourceLines)
% COUNT_UNCLEAR_SINGLE_LETTERS Count unclear one-letter assignments.
%
% Syntax:
%   unclearSingleLetterCount = count_unclear_single_letters(sourceLines)
%
% Description:
%   Excludes approved matrix names, i/j/k loop indices, and Cartesian
%   coordinate or velocity components x/y/z and u/v/w.
%
% Inputs:
%   sourceLines - MATLAB source lines, string array [-].
%
% Outputs:
%   unclearSingleLetterCount - Number of unclear assignments [-].
%
% Governing Equations / Theory:
%   Uses a deterministic case-sensitive regular-expression match.
%
% References:
%   - CRESTU engineering code standard.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

    %% Stage 1: Match unapproved one-letter assignments

    expression = '^\s*(for\s+)?[bcefghlnoqrsp]\s*=';
    unclearSingleLetterCount = sum(~cellfun('isempty', ...
        regexp(cellstr(sourceLines), expression, 'once')));
end

function compactControlCount = count_compact_controls(sourceLines)
% COUNT_COMPACT_CONTROLS Count control statements with inline execution.
%
% Syntax:
%   compactControlCount = count_compact_controls(sourceLines)
%
% Description:
%   Counts lines that place a control header and executable body together.
%
% Inputs:
%   sourceLines - MATLAB source lines, string array [-].
%
% Outputs:
%   compactControlCount - Number of compact control lines [-].
%
% Governing Equations / Theory:
%   Uses a deterministic regular-expression match.
%
% References:
%   - CRESTU engineering code standard.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

    %% Stage 1: Match compact control syntax

    expression = ['^\s*((if|elseif|for|while)\b.*,' ...
        '\s*[A-Za-z_]\w*\s*=|case\b.*,\s*[A-Za-z_]\w*|otherwise\s*,)'];
    compactControlCount = sum(~cellfun('isempty', ...
        regexp(cellstr(sourceLines), expression, 'once')));
end

function nonAsciiLineCount = count_non_ascii_lines(sourceLines)
% COUNT_NON_ASCII_LINES Count source lines that contain non-ASCII text.
%
% Syntax:
%   nonAsciiLineCount = count_non_ascii_lines(sourceLines)
%
% Description:
%   Detects characters outside the ASCII range used by project source files.
%
% Inputs:
%   sourceLines - MATLAB source lines, string array [-].
%
% Outputs:
%   nonAsciiLineCount - Number of source lines with non-ASCII text [-].
%
% Governing Equations / Theory:
%   Uses a deterministic regular-expression match.
%
% References:
%   - CRESTU engineering code standard.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

    %% Stage 1: Match characters outside the ASCII range

    nonAsciiLineCount = sum(~cellfun('isempty', ...
        regexp(cellstr(sourceLines), '[^\x00-\x7F]', 'once')));
end
