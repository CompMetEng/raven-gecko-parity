d = getenv('P5_DIR'); res = struct();
orig = pr_common('model');
cases = {'met','rxn','gene','space'};
for i = 1:numel(cases)
    f = fullfile(d, ['empty_' cases{i} '.yml']); r = struct('file', f);
    txt = strsplit(fileread(f), newline);
    if strcmp(cases{i}, 'gene')
        r.geneNameKeyLinesInFile = sum(contains(txt, 'name: ""')) ;
    end
    try
        evalc('m = readYAMLmodel(f);');
        r.loaded = true;
        r.metNames1 = m.metNames(1); r.rxnNames1 = m.rxnNames(1);
        if isfield(m, 'geneShortNames'), r.geneShortNames1 = m.geneShortNames(1); end
        k = 1:5;
        r.first5.mets = m.mets(k); r.first5.metNames = m.metNames(k); r.first5.metFormulas = m.metFormulas(k);
        r.first5Orig.mets = orig.mets(k); r.first5Orig.metNames = orig.metNames(k); r.first5Orig.metFormulas = orig.metFormulas(k);
        r.metsEqualOrig = isequal(m.mets, orig.mets); r.metFormulasEqualOrig = isequal(m.metFormulas, orig.metFormulas);
        r.nonEmptyNameCount = [sum(~cellfun(@isempty, m.metNames)) sum(~cellfun(@isempty, orig.metNames))];
        try
            msg = evalc('checkModelStruct(m, false);');
            r.checkModelStruct = msg;
        catch e2, r.checkModelStruct = struct('error', pr_common('err', e2)); end
    catch e
        r.loaded = false; r.error = pr_common('err', e);
    end
    res.(cases{i}) = r;
end
% single-space name written by RAVEN
m = orig; m.metNames{1} = ' ';
f = [tempname '.yml'];
try
    evalc('writeYAMLmodel(m, f);');
    t = strsplit(fileread(f), newline); i1 = find(contains(t, 'metabolites:'), 1);
    res.writeSingleSpace.linesAfterMetabolites = t(i1:min(i1+4, numel(t)));
catch e, res.writeSingleSpace.error = pr_common('err', e); end
pr_common('save', '5_emptyNames', res);
