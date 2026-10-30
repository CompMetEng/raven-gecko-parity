function results = ec_model_load_light_humangem(ctx)
% MATLAB side of the light Human-GEM ecModel loading scenario; mirrors run.py. loadEcModel
% reads the light ecModel of GECKO's light_ecModel tutorial as published and as repaired, and
% the repaired model is summarised field by field. See scenario.yml for what is compared:
% lists keep the order of the model, the stoichiometry and the coupling matrix are sorted
% multi-key records joined on char(1), an infinite bound or a NaN is a class and a zero.
%
% loadEcModel joins its filename onto <adapter path>/models/. A copy of the adapter with its
% path pointed at a temporary folder routes both files there; ModelAdapter is a value class.

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
params = adapter.getParameters();
workdir = tempname;
modelsDir = fullfile(workdir, 'models');
mkdir(modelsDir);
cleanup = onCleanup(@() rmdir(workdir, 's'));
ioAdapter = adapter;
ioAdapter.params.path = workdir;

text = fileread(fullfile(params.path, 'models', char(ctx.inputs.model_file)));
write_text(fullfile(modelsDir, 'published.yml'), text);
% The file has 52 ec reactions written as an `enzymes: !!omap` entry with no content (the
% reactions that carry the standard kcat). as_published records what each reader does with
% the file as it is; `repaired` is the summary of a copy without those empty entries, which
% both readers can load.
repaired = regexprep(text, '(?m)^ *- enzymes: !!omap\n(?= {0,4}\S)', '');
write_text(fullfile(modelsDir, 'repaired.yml'), repaired);

warnState = warning('off', 'all');
restoreWarnings = onCleanup(@() warning(warnState));

results.as_published = attempt(@() counts_of('published.yml', ioAdapter));
results.repaired = summary(loadEcModel('repaired.yml', ioAdapter));

end


function write_text(path, text)
fid = fopen(path, 'w', 'n', 'UTF-8');
fprintf(fid, '%s', text);
fclose(fid);
end


function out = attempt(fn)
% A raise is a result on both sides.
try
    out = struct('raised', false, 'value', fn());
catch err
    fprintf('case raised: %s\n', err.message);
    out = struct('raised', true, 'value', struct());
end
end


function out = counts_of(name, ioAdapter)
model = loadEcModel(name, ioAdapter);
out = struct('reactions', numel(model.rxns), 'metabolites', numel(model.mets), ...
    'genes', numel(model.genes), 'ec_rxns', numel(model.ec.rxns), 'enzymes', numel(model.ec.enzymes));
end


function out = summary(model)
out.counts = struct( ...
    'reactions', numel(model.rxns), ...
    'metabolites', numel(model.mets), ...
    'genes', numel(model.genes), ...
    'ec_rxns', numel(model.ec.rxns), ...
    'enzymes', numel(model.ec.enzymes));
out.model_id = char(model.id);
out.gecko_light = logical(model.ec.geckoLight);
out.reaction_order = reshape(model.rxns, 1, []);
out.metabolite_order = reshape(model.mets, 1, []);
out.gene_order = reshape(model.genes, 1, []);
out.compartments = reshape(sort(model.comps), 1, []);
out.reactions = reactions(model);
out.gene_rules = gene_rules(model);
out.metabolites = metabolites(model);
out.stoichiometry = pairs(model.S.', model.rxns, model.mets, {'reaction', 'metabolite'}, true);
out.ec = ec_data(model.ec);
end


function out = reactions(model)
n = numel(model.rxns);
out = cell(1, n);
for i = 1:n
    [lowerValue, lowerKind] = bound(model.lb(i));
    [upperValue, upperKind] = bound(model.ub(i));
    out{i} = struct( ...
        'id',                    model.rxns{i}, ...
        'name',                  char(model.rxnNames{i}), ...
        'lower_bound',           lowerValue, ...
        'lower_kind',            lowerKind, ...
        'upper_bound',           upperValue, ...
        'upper_kind',            upperKind, ...
        'objective_coefficient', double(model.c(i)), ...
        'subsystem',             {names(subsystem_of(model, i))}, ...
        'ec_codes',              ec_codes_of(model, i), ...
        'genes',                 {names(model.genes(model.rxnGeneMat(i, :) ~= 0))});
end
end


function out = gene_rules(model)
n = numel(model.rxns);
out = cell(1, n);
for i = 1:n
    rule = '';
    if isfield(model, 'grRules') && numel(model.grRules) >= i
        rule = char(model.grRules{i});
    end
    out{i} = struct('id', model.rxns{i}, 'rule', rule);
end
end


function out = metabolites(model)
n = numel(model.mets);
out = cell(1, n);
for i = 1:n
    charge = NaN;
    if isfield(model, 'metCharges') && numel(model.metCharges) >= i
        charge = double(model.metCharges(i));
    end
    formula = '';
    if isfield(model, 'metFormulas') && numel(model.metFormulas) >= i
        formula = char(model.metFormulas{i});
    end
    out{i} = struct( ...
        'id',          model.mets{i}, ...
        'name',        char(model.metNames{i}), ...
        'compartment', model.comps{model.metComps(i)}, ...
        'formula',     formula, ...
        'charge_set',  ~isnan(charge), ...
        'charge',      zero_for_nan(charge));
end
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


function out = ec_data(ec)
out.rxns     = reshape(ec.rxns(:), 1, []);
out.genes    = reshape(ec.genes(:), 1, []);
out.enzymes  = reshape(ec.enzymes(:), 1, []);
[mw, mwKind] = numbers(ec.mw);
out.mw       = mw;
out.mw_kind  = mwKind;
out.sequence = reshape(blanks_for_empty(ec.sequence), 1, []);
[concs, concsKind] = numbers(ec.concs);
out.concs      = concs;
out.concs_kind = concsKind;
out.eccodes  = reshape(blanks_for_empty(ec.eccodes), 1, []);
out.kcat     = num2cell(reshape(double(ec.kcat(:)), 1, []));
out.source   = reshape(blanks_for_empty(ec.source), 1, []);
out.notes    = reshape(blanks_for_empty(ec.notes), 1, []);
out.coupling = pairs(ec.rxnEnzMat, ec.rxns, ec.enzymes, {'reaction', 'enzyme'}, true);
end


function [values, kinds] = numbers(x)
% A number vector as a cell of values and a cell of kinds: NaN and infinity are classes.
x = double(x(:)).';
values = cell(1, numel(x));
kinds = cell(1, numel(x));
for k = 1:numel(x)
    if isnan(x(k))
        values{k} = 0;
        kinds{k} = 'nan';
    else
        [values{k}, kinds{k}] = bound(x(k));
    end
end
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
