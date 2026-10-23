function results = complex_data_edge_cases_ectestgem(ctx)
% MATLAB side of the complex-data edge-case scenario.
%
% Must return the shape run.py returns. Conventions from the earlier gecko scenarios:
% rxnEnzMat is emitted as (reaction, enzyme, coefficient) records sorted on the keys joined
% with char(1); NaN is never emitted (a flag and a zero); an error is a flag, not a message;
% a missing data file fails the run instead of being caught.

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);

% --- complexes
names = sort(cellstr(ctx.inputs.complex_cases));
complexes = struct();
for k = 1:numel(names)
    if strcmp(names{k}, 'fixture')
        file = fullfile(ctx.inputs.fixture_dir, 'data', 'ComplexPortal.json');
    else
        file = fullfile(ctx.inputs.complex_dir, [names{k} '.json']);
    end
    if ~isfile(file)
        error('complex_data_edge_cases_ectestgem:noFile', 'missing %s', file);
    end
    perFlavour = struct();
    flavours = {'full', false; 'light', true};
    for f = 1:size(flavours, 1)
        ec = makeEcModel(model, flavours{f, 2}, adapter);
        try
            ec = applyComplexData(ec, file, adapter, false);
            perFlavour.(flavours{f, 1}) = struct('raised', false, ...
                'coupling', {pairs(ec.ec.rxnEnzMat, ec.ec.rxns, ec.ec.enzymes)});
        catch
            perFlavour.(flavours{f, 1}) = struct('raised', true, 'coupling', {{}});
        end
    end
    complexes.(names{k}) = perFlavour;
end

% --- f factor
ecForF = makeEcModel(model, false, adapter);
fCases = ctx.inputs.f_factor_cases;
names = sort(fieldnames(fCases));
fFactor = struct();
for k = 1:numel(names)
    c = fCases.(names{k});
    protData = struct('uniprotIDs', {cellstr(c.ids)}, 'abundances', double(c.abundances));
    args = {'protData', protData, 'modelAdapter', adapter};
    if c.empty_enzymes
        args = [args, {'enzymes', {}}]; %#ok<AGROW>
    end
    try
        value = double(calculateFfactor(ecForF, args{:}));
        isNan = isnan(value);
        if isNan
            value = 0;
        end
        fFactor.(names{k}) = struct('raised', false, 'is_nan', isNan, 'value', value);
    catch
        fFactor.(names{k}) = struct('raised', true, 'is_nan', false, 'value', 0);
    end
end

% --- fill concs
fillCases = ctx.inputs.fill_cases;
names = sort(fieldnames(fillCases));
fillConcs = struct();
for k = 1:numel(names)
    c = fillCases.(names{k});
    protData = struct('uniprotIDs', {cellstr(c.ids)}, 'abundances', double(c.abundances));
    ec = makeEcModel(model, false, adapter);
    try
        ec = fillEnzConcs(ec, protData, double(c.data_col));
        fillConcs.(names{k}) = struct('raised', false, 'concs', {concs(ec)});
    catch
        fillConcs.(names{k}) = struct('raised', true, 'concs', {{}});
    end
end

results.complexes = complexes;
results.f_factor = fFactor;
results.fill_concs = fillConcs;

end


function out = concs(model)
% ec.concs paired with the enzyme it belongs to, sorted by enzyme id.
[enzymeIds, order] = sort(model.ec.enzymes(:));
values = double(model.ec.concs(:));
values = values(order);

out = cell(1, numel(enzymeIds));
for k = 1:numel(enzymeIds)
    hasConc = ~isnan(values(k));
    if hasConc
        conc = values(k);
    else
        conc = 0;
    end
    out{k} = struct('enzyme', enzymeIds{k}, 'has_conc', hasConc, 'conc', conc);
end
end


function out = pairs(matrix, rowIds, colIds)
% Non-zero entries of ec.rxnEnzMat as a sorted list of records (char(1) keeps Python's
% tuple order).
[rows, cols, values] = find(matrix);
keys = cell(numel(rows), 1);
for k = 1:numel(rows)
    keys{k} = [rowIds{rows(k)} char(1) colIds{cols(k)}];
end
[~, order] = sort(keys);

out = cell(1, numel(order));
for k = 1:numel(order)
    i = order(k);
    out{k} = struct( ...
        'reaction',    rowIds{rows(i)}, ...
        'enzyme',      colIds{cols(i)}, ...
        'coefficient', double(values(i)));
end
end
