function hashText = sha256_file_hash(filename)
% SHA256_FILE_HASH Compute a SHA-256 digest from the exact file bytes.
%
% Syntax:
%   hashText = sha256_file_hash(filename)
%
% Inputs:
%   filename - Existing file path [-].
%
% Outputs:
%   hashText - Lowercase hexadecimal SHA-256 digest [64 characters].

    arguments
        filename {mustBeTextScalar}
    end

    %% 阶段 1: 验证并读取文件

    assert(isfile(filename), 'CRESTU:HashFileMissing', ...
        'File selected for hashing does not exist: %s', filename);
    fileIdentifier = fopen(filename, 'rb');
    assert(fileIdentifier >= 0, 'CRESTU:HashFileOpen', ...
        'Unable to open file for hashing: %s', filename);
    fileCleanup = onCleanup(@() fclose(fileIdentifier));
    fileBytes = fread(fileIdentifier, Inf, '*uint8');
    clear fileCleanup

    %% 阶段 2: 计算文件摘要

    hashText = sha256_hash(fileBytes);
end
