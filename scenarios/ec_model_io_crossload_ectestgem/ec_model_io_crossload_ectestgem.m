function results = ec_model_io_crossload_ectestgem(ctx)
% MATLAB side of the ecModel cross-load scenario.
%
% loadEcModel reads the file saveEcModel wrote and the file save_ec_model wrote (both
% vendored under data/), and each loaded model is summarised the way run.py does. See
% scenario.yml for why the files are fixtures and how they are kept current. The
% conventions of ec_model_io_ectestgem carry over: an infinite bound is a class and a zero;
% ec.rxns / ec.genes / ec.enzymes stay in their own order; stoichiometry and the coupling
% matrix are sorted multi-key records joined on char(1); ec.concs is left out.

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
dataDir = fullfile(ctx.inputs.scenario_root, 'data');
model = loadConventionalGEM('modelAdapter', adapter);
ecModel = makeEcModel(model, false, adapter);
% R5 becomes an ec reaction without enzymes, as a standard kcat row is.
row = find(strcmp(ecModel.ec.rxns, 'R5'));
ecModel.ec.rxnEnzMat(row, :) = 0;
ecModel.ec.source{row} = 'standard';

workdir = tempname;
mkdir(fullfile(workdir, 'models'));
cleanup = onCleanup(@() rmdir(workdir, 's'));

% saveEcModel and loadEcModel join their filename onto <adapter path>/models/. A throwaway
% copy of the adapter with its path pointed at a temp directory routes both there instead
% of into the fixture's own models/ folder. ModelAdapter is a value class.
ioAdapter = adapter;
ioAdapter.params.path = workdir;
modelsDir = fullfile(workdir, 'models');

saveEcModel(ecModel, 'own.yml', ioAdapter);
if ~isequal(without_date(fullfile(modelsDir, 'own.yml')), ...
        without_date(fullfile(dataDir, 'ecModel.matlab.yml')))
    error('saveEcModel no longer writes data/ecModel.matlab.yml; regenerate the fixtures with regenerate_fixtures.py');
end

copyfile(fullfile(dataDir, 'ecModel.python.yml'), fullfile(modelsDir, 'from_python.yml'));
copyfile(fullfile(dataDir, 'ecModel.matlab.yml'), fullfile(modelsDir, 'from_matlab.yml'));

results.file_written_by_python = load_summary('from_python.yml', ioAdapter);
results.file_written_by_matlab = load_summary('from_matlab.yml', ioAdapter);

end


function lines = without_date(path)
text = strrep(fileread(path), sprintf('\r\n'), sprintf('\n'));
lines = strsplit(text, sprintf('\n'), 'CollapseDelimiters', false);
lines = lines(~startsWith(lines, '  - date:'));
end


function out = load_summary(filename, adapter)
% A loader that refuses the file is the result, not a failure of the scenario.
try
    out = summary(loadEcModel(filename, adapter));
catch err
    out = struct('loaded', false, 'error', err.identifier);
    if isempty(err.identifier)
        out.error = 'error';
    end
end
end


function out = summary(model)
out.loaded = true;
out.model_id = model.id;
out.description = '';
out.geckopy_version = '';
if isfield(model, 'annotation')
    if isfield(model.annotation, 'note')
        out.description = char(model.annotation.note);
    end
    if isfield(model.annotation, 'geckopy_version')
        out.geckopy_version = char(model.annotation.geckopy_version);
    end
end
out.n_reactions = numel(model.rxns);
out.n_metabolites = numel(model.mets);
out.n_genes = numel(model.genes);

[rxnIds, rxnOrder] = sort(model.rxns(:));
reactions = cell(1, numel(rxnIds));
for k = 1:numel(rxnIds)
    i = rxnOrder(k);
    [lowerValue, lowerKind] = bound(model.lb(i));
    [upperValue, upperKind] = bound(model.ub(i));
    reactions{k} = struct( ...
        'id',                     rxnIds{k}, ...
        'lower_bound',            lowerValue, ...
        'lower_kind',             lowerKind, ...
        'upper_bound',            upperValue, ...
        'upper_kind',             upperKind, ...
        'objective_coefficient',  double(model.c(i)), ...
        'subsystem',              {names(subsystem_of(model, i))}, ...
        'ec_codes',               ec_codes_of(model, i), ...
        'sbo',                    sbo_of(model, 'rxnMiriams', i));
end
out.reactions = reactions;

[metIds, metOrder] = sort(model.mets(:));
metabolites = cell(1, numel(metIds));
for k = 1:numel(metIds)
    i = metOrder(k);
    charge = NaN;
    if isfield(model, 'metCharges') && numel(model.metCharges) >= i
        charge = double(model.metCharges(i));
    end
    metabolites{k} = struct( ...
        'id',          metIds{k}, ...
        'name',        model.metNames{i}, ...
        'compartment', model.comps{model.metComps(i)}, ...
        'charge_set',  ~isnan(charge), ...
        'charge',      zero_for_nan(charge), ...
        'sbo',         sbo_of(model, 'metMiriams', i));
end
out.metabolites = metabolites;

[geneIds, geneOrder] = sort(model.genes(:));
genes = cell(1, numel(geneIds));
for k = 1:numel(geneIds)
    i = geneOrder(k);
    name = '';
    if isfield(model, 'geneShortNames') && numel(model.geneShortNames) >= i
        name = char(model.geneShortNames{i});
    end
    genes{k} = struct('id', geneIds{k}, 'name', name, 'sbo', sbo_of(model, 'geneMiriams', i));
end
out.genes = genes;

out.gene_associations = pairs(model.rxnGeneMat, model.rxns, model.genes, ...
    {'reaction', 'gene'}, false);
out.stoichiometry = pairs(model.S.', model.rxns, model.mets, ...
    {'reaction', 'metabolite'}, true);

out.gecko_light = logical(model.ec.geckoLight);
out.ec = ec_data(model.ec);
end


function value = zero_for_nan(x)
value = x;
if isnan(x)
    value = 0;
end
end


function out = names(values)
% A sorted row of the non-empty strings in values, which may be a char, a cell or empty.
if isempty(values)
    out = {};
    return
end
if ischar(values) || isstring(values)
    values = cellstr(values);
end
values = values(:).';
values = values(~cellfun(@isempty, values));
out = sort(values);
end


function value = subsystem_of(model, i)
value = {};
if isfield(model, 'subSystems') && numel(model.subSystems) >= i
    value = model.subSystems{i};
end
end


function value = ec_codes_of(model, i)
value = '';
if isfield(model, 'eccodes') && numel(model.eccodes) >= i && ~isempty(model.eccodes{i})
    value = char(model.eccodes{i});
end
end


function value = sbo_of(model, field, i)
% The SBO terms of one entity, ';'-joined, read from its MIRIAM annotation.
value = '';
if ~isfield(model, field) || numel(model.(field)) < i || isempty(model.(field){i})
    return
end
miriam = model.(field){i};
if isstruct(miriam) && isfield(miriam, 'name')
    hit = strcmp(miriam.name, 'sbo');
    found = names(miriam.value(hit));
    value = strjoin(found, ';');
end
end


function out = ec_data(ec)
out.rxns     = reshape(ec.rxns(:), 1, []);
out.genes    = reshape(ec.genes(:), 1, []);
out.enzymes  = reshape(ec.enzymes(:), 1, []);
out.mw       = reshape(double(ec.mw(:)), 1, []);
out.sequence = reshape(ec.sequence(:), 1, []);
out.eccodes  = reshape(blanks_for_empty(ec.eccodes), 1, []);
out.kcat     = reshape(double(ec.kcat(:)), 1, []);
out.source   = reshape(blanks_for_empty(ec.source), 1, []);
out.coupling = pairs(ec.rxnEnzMat, ec.rxns, ec.enzymes, {'reaction', 'enzyme'}, true);
end


function [value, kind] = bound(x)
if isinf(x) && x > 0
    value = 0;
    kind = '+inf';
elseif isinf(x) && x < 0
    value = 0;
    kind = '-inf';
else
    value = double(x);
    kind = 'finite';
end
end


function out = pairs(matrix, rowIds, colIds, names, withCoefficient)
[rows, cols, values] = find(matrix);
keys = cell(numel(rows), 1);
for k = 1:numel(rows)
    keys{k} = [rowIds{rows(k)} char(1) colIds{cols(k)}];
end
[~, order] = sort(keys);

out = cell(1, numel(order));
for k = 1:numel(order)
    i = order(k);
    record = struct(names{1}, rowIds{rows(i)}, names{2}, colIds{cols(i)});
    if withCoefficient
        record.coefficient = double(values(i));
    end
    out{k} = record;
end
end


function out = blanks_for_empty(values)
out = values(:).';
for k = 1:numel(out)
    if isempty(out{k})
        out{k} = '';
    end
end
end
