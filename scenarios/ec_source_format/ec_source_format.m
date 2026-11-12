function results = ec_source_format(ctx)
% MATLAB side of the ec.source format scenario; mirrors run.py. parseKcatSource and
% formatKcatSource on the strings and calls of inputs, and the source that setKcatForReactions
% writes on ecTestGEM. A case that raises is a result, not a failure.

inputs = ctx.inputs;
results = struct();

texts = inputs.parse;
for i = 1:numel(texts)
    results.(sprintf('parse_%02d', i - 1)) = run_case(@() parse_record(texts{i}));
end

specs = inputs.format;
for i = 1:numel(specs)
    results.(sprintf('format_%02d', i - 1)) = run_case(@() format_record(specs{i}));
end

results.set_kcat_manual = run_case(@() set_kcat_record(ctx));

end


function out = parse_record(text)
info = parseKcatSource(text);
keys = sort(fieldnames(info.detail));
detail = cell(1, numel(keys));
for i = 1:numel(keys)
    value = info.detail.(keys{i});
    if ischar(value)
        detail{i} = struct('key', keys{i}, 'kind', 'text', 'number', num(NaN), 'text', value);
    else
        detail{i} = struct('key', keys{i}, 'kind', 'number', 'number', num(value), 'text', '');
    end
end
out = struct('origin', info.origin, 'detail', {detail}, 'tuned', {reshape(info.tuned, 1, [])}, ...
    'prior', num(info.prior));
end


function out = format_record(spec)
origin = spec.origin;
detail = [];
if isfield(spec, 'detail') && ~isempty(spec.detail)
    detail = spec.detail;
end
tuned = {};
if isfield(spec, 'tuned') && ~isempty(spec.tuned)
    tuned = cellstr(spec.tuned);
end
prior = NaN;
if isfield(spec, 'prior') && ~isempty(spec.prior)
    prior = spec.prior;
end
out = struct('text', formatKcatSource(origin, detail, tuned, prior));
end


function out = set_kcat_record(ctx)
inputs = ctx.inputs;
spec = inputs.set_kcat;
adapter = ModelAdapterManager.getAdapter(inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);
index = find(strcmp(ecModel.ec.rxns, spec.rxn));
ecModel.ec.source{index} = spec.source_before;
others = ecModel.ec.source;
others(index) = [];
ecModel = setKcatForReactions(ecModel, spec.rxn, double(spec.kcat));
after = ecModel.ec.source;
after(index) = [];
out = struct('source', ecModel.ec.source{index}, 'kcat', num(ecModel.ec.kcat(index)), ...
    'others_unchanged', isequal(others, after));
end


% --------------------------------------------------------------------------- %
% shared helpers (same in every probe scenario of this family)
% --------------------------------------------------------------------------- %

function out = num(x)
% A number as a record, so that NaN is a value both sides can state.
x = double(x);
if isempty(x)
    x = NaN;
end
out = struct('is_nan', logical(isnan(x)), 'value', double(nan2zero(x)));
end


function x = nan2zero(x)
x(isnan(x)) = 0;
end


function out = run_case(fn)
% Result of one case; a raise is a result on both sides, not a failure.
try
    value = fn();
catch
    out = struct('raised', true, 'value', struct());
    return
end
out = struct('raised', false, 'value', value);
end
