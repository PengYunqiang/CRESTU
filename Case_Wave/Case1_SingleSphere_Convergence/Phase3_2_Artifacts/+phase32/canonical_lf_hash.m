function hashText = canonical_lf_hash(filename)
% CANONICAL_LF_HASH Hash text after in-memory CRLF/CR to LF conversion.

    assert(isfile(filename), 'CRESTU:Phase32CanonicalFileMissing', ...
        'Cannot hash a missing file: %s', filename);
    fileID = fopen(filename, 'rb');
    assert(fileID >= 0, 'CRESTU:Phase32CanonicalFileOpen', ...
        'Cannot open file for hashing: %s', filename);
    cleanup = onCleanup(@() fclose(fileID));
    bytes = fread(fileID, Inf, '*uint8')';
    clear cleanup
    textValue = native2unicode(bytes, 'UTF-8');
    textValue = strrep(textValue, [char(13), newline], newline);
    textValue = strrep(textValue, char(13), newline);
    hashText = sha256_hash(uint8(unicode2native(textValue, 'UTF-8')));
end
