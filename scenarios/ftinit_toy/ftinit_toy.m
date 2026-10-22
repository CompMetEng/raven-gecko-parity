function results = ftinit_toy(ctx)
% MATLAB side of the prepINITModel and ftINIT scenario on RAVEN's ten-reaction test model.
%
% Each case prepares the model (prepINITModel) and extracts a model from the expression levels of one
% tissue (ftINIT with 'transcrData'). Reaction and gene lists are sorted so that only the content is
% compared.

inputs = ctx.inputs;
previousSolver = getpref('RAVEN', 'solver', '');
restoreSolver = onCleanup(@() restore_solver(previousSolver));
setRavenSolver(char(inputs.matlab_solver));

cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    results.(char(cases{k}.name)) = run_case(inputs, cases{k});
end
end


function out = run_case(inputs, c)
model = build_model(inputs, logical(c.overrides));
if logical(c.with_task)
    taskStruct = make_task(inputs.task);
else
    taskStruct = [];
end
spont = as_cellstr(c.spontaneous);
evalc('prepData = prepINITModel(model, taskStruct, ''spontRxnNames'', spont, ''extComp'', char(inputs.ext_comp));');

genes = sort(model.genes(:));
levels = zeros(numel(genes), 1);
for i = 1:numel(genes)
    levels(i) = double(inputs.levels.(genes{i}));
end
arrayData = struct('genes', {genes}, 'tissues', {{char(inputs.tissue)}}, 'celltypes', {{''}}, ...
    'levels', levels, 'threshold', double(inputs.threshold));

if strcmp(char(c.series), 'full')
    steps = getINITSteps([], 'full');
else
    steps = getINITSteps();
end
evalc(['res = ftINIT(prepData, char(inputs.tissue), [], [], ''transcrData'', arrayData, ' ...
    '''INITSteps'', steps, ''removeGenes'', logical(c.remove_genes));']);

ref = prepData.refModel;
masks = struct( ...
    'exchange', prepData.toIgnoreExch, 'import', prepData.toIgnoreImportRxns, ...
    'simple_transport', prepData.toIgnoreSimpleTransp, 'advanced_transport', prepData.toIgnoreAdvTransp, ...
    'spontaneous', prepData.toIgnoreSpont, 'no_gpr', prepData.toIgnoreAllWithoutGPRs);
prep = struct('ref_reactions', {sort_cellstr(ref.rxns)}, ...
    'n_min_reactions', numel(prepData.minModel.rxns), ...
    'essential_reactions', {sort_cellstr(prepData.essentialRxns)});
for f = fieldnames(masks)'
    mask = masks.(f{1});
    if numel(mask) == numel(ref.rxns)
        prep.(f{1}) = sort_cellstr(ref.rxns(logical(mask)));
    else
        prep.(f{1}) = sprintf('mask of %d elements for %d reactions', numel(mask), numel(ref.rxns));
    end
end

objective = char(inputs.objective);
idx = find(strcmp(res.rxns, objective));
flux = 0;
if ~isempty(idx)
    m = res;
    m.c = zeros(numel(m.rxns), 1);
    m.c(idx) = 1;
    sol = solveLP(m);
    if ~isempty(sol.x)
        flux = sol.x(idx);
    else
        flux = NaN;
    end
end

[~, order] = sort(res.rxns);
rules = cell(1, numel(order));
for j = 1:numel(order)
    rule = res.grRules{order(j)};
    clauses = {};
    if ~isempty(rule)
        dnf = grRuleToDNF(rule);
        clauses = cell(1, numel(dnf));
        for q = 1:numel(dnf)
            clauses{q} = sort(dnf{q}(:))';
        end
        [~, o] = sort(cellfun(@(x) strjoin(x, ','), clauses, 'UniformOutput', false));
        clauses = clauses(o);
    end
    rules{j} = struct('reaction', res.rxns{order(j)}, 'clauses', {clauses});
end
[~, bo] = sort(res.rxns);
bounds = cell(1, numel(bo));
for j = 1:numel(bo)
    bounds{j} = struct('reaction', res.rxns{bo(j)}, 'lower_bound', double(res.lb(bo(j))), 'upper_bound', double(res.ub(bo(j))));
end
out = struct('prep', prep, 'model', struct( ...
    'bounds', {bounds}, 'reactions', {sort_cellstr(res.rxns)}, 'genes', {sort_cellstr(res.genes)}, 'rules', {rules}, ...
    'has_objective', ~isempty(idx), 'objective_flux', flux));
end


function model = build_model(inputs, overrides)
mets = as_list(inputs.metabolites);
rxns = as_list(inputs.reactions);
rules = cell(numel(rxns), 1);
for i = 1:numel(rxns)
    rules{i} = char(rxns{i}.rule);
    if overrides && isfield(inputs.rule_overrides, rxns{i}.id)
        rules{i} = char(inputs.rule_overrides.(rxns{i}.id));
    end
end
genes = {};
for i = 1:numel(rules)
    genes = [genes, regexp(rules{i}, 'G\d+', 'match')]; %#ok<AGROW>
end
genes = unique(genes, 'stable')';

model = struct();
model.id = 'testModel';
model.rxns = {};
model.S = [];
model.rev = [];
model.mets = cellfun(@(m) char(m.id), mets, 'UniformOutput', false)';
model.metNames = cellfun(@(m) char(m.name), mets, 'UniformOutput', false)';
metComps = cellfun(@(m) char(m.compartment), mets, 'UniformOutput', false)';
model.comps = unique(metComps, 'stable');
model.compNames = model.comps;
[~, model.metComps] = ismember(metComps, model.comps);
model.genes = genes;
model.grRules = {};
model.rxnGeneMat = [];

toAdd = struct();
toAdd.rxns = cellfun(@(r) char(r.id), rxns, 'UniformOutput', false)';
toAdd.equations = cellfun(@(r) char(r.equation), rxns, 'UniformOutput', false)';
toAdd.grRules = rules;
evalc('model = addRxns(model, toAdd, 3);');
model.lb = cellfun(@(r) double(r.lb), rxns)';
model.ub = cellfun(@(r) double(r.ub), rxns)';
model.c = double(strcmp(model.rxns, char(inputs.objective)));
model.rxnNames = model.rxns;
model.b = zeros(numel(model.mets), 1);
end


function task = make_task(spec)
task = struct();
task.id = char(spec.id);
task.description = char(spec.description);
task.shouldFail = false;
task.printFluxes = false;
task.comments = '';
ins = as_list(spec.inputs);
outs = as_list(spec.outputs);
task.inputs = cellfun(@(i) char(i.metabolite), ins, 'UniformOutput', false)';
task.LBin = cellfun(@(i) double(i.lb), ins)';
task.UBin = inf(numel(ins), 1);
task.outputs = cellfun(@(o) char(o.metabolite), outs, 'UniformOutput', false)';
task.LBout = cellfun(@(o) double(o.lb), outs)';
task.UBout = cellfun(@(o) double(o.ub), outs)';
task.equations = {};
task.LBequ = [];
task.UBequ = [];
task.changed = {};
task.LBrxn = {};
task.UBrxn = {};
end


function restore_solver(previous)
if ~isempty(previous)
    setRavenSolver(previous);
end
end


function out = sort_cellstr(values)
if isempty(values)
    out = {};
    return
end
out = reshape(sort(values(:)), 1, []);
end


function out = as_cellstr(value)
if isempty(value)
    out = {};
elseif ischar(value)
    out = {value};
else
    out = reshape(cellfun(@char, value, 'UniformOutput', false), 1, []);
end
end


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end
