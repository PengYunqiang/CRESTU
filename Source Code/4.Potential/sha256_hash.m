function hashText = sha256_hash(inputData)
% SHA256_HASH Compute a lowercase SHA-256 digest for bytes or UTF-8 text.
%
% Syntax:
%   hashText = sha256_hash(inputData)
%
% Inputs:
%   inputData - Byte vector, character vector, or string scalar [-].
%
% Outputs:
%   hashText  - Lowercase hexadecimal SHA-256 digest [64 characters].

    %% 阶段 1: 规范化待哈希数据

    if ischar(inputData) || (isstring(inputData) && isscalar(inputData))
        inputBytes = unicode2native(char(inputData), 'UTF-8');
    elseif isa(inputData, 'uint8')
        inputBytes = inputData;
    else
        error('CRESTU:HashInputType', ...
            'SHA-256 input must be uint8 bytes or scalar text.');
    end

    inputBytes = reshape(inputBytes, 1, []);

    %% 阶段 2: 计算 SHA-256 摘要

    digestEngine = java.security.MessageDigest.getInstance('SHA-256');
    digestEngine.update(inputBytes);
    signedDigest = digestEngine.digest();
    digestBytes = typecast(signedDigest, 'uint8');
    hashText = lower(reshape(dec2hex(digestBytes, 2).', 1, []));
end
