function results = reporter_metabolites_toy(ctx)
% MATLAB side of the reporter-metabolites scenario on a hand-written model.
%
% The result of each case is the list of tests the function returns, each with the metabolites in the
% order returned and their scores. 'only up' and 'only down' are reported as 'up' and 'down'.

inputs = ctx.inputs;
reactions = as_list(inputs.reactions);
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    model = toy_model(reactions);
    genes = as_cellstr(c.genes);
    pvalues = double(c.pvalues(:));
    args = {};
    if ~isempty(c.fold_changes)
        args = {'geneFoldChanges', double(c.fold_changes(:))};
    end
    rep = reporterMetabolites(model, genes(:), pvalues, args{:});
    tests = cell(1, numel(rep));
    for t = 1:numel(rep)
        r = rep(t);
        name = strrep(strrep(r.test, 'only ', ''), ' ', '_');
        tests{t} = struct('test', name, ...
            'mets', {reshape(r.mets, 1, [])}, ...
            'z_scores', {reshape(double(r.metZScores), 1, [])}, ...
            'p_values', {reshape(double(r.metPValues), 1, [])}, ...
            'n_genes', {reshape(double(r.metNGenes), 1, [])}, ...
            'mean_z', {reshape(double(r.meanZ), 1, [])}, ...
            'std_z', {reshape(double(r.stdZ), 1, [])}, ...
            'std_is_nan', {reshape(isnan(r.stdZ), 1, [])});
    end
    results.(char(c.name)) = struct('n_results', numel(rep), 'tests', {tests});
end
end


function model = toy_model(reactions)
n = numel(reactions);
mets = {};
genes = {};
for i = 1:n
    mets = [mets, as_cellstr(reactions{i}.substrates), as_cellstr(reactions{i}.products)]; %#ok<AGROW>
    genes = [genes, rule_genes(char(reactions{i}.rule))]; %#ok<AGROW>
end
mets = unique(mets, 'stable');
genes = unique(genes, 'stable');
S = zeros(numel(mets), n);
rxnGeneMat = zeros(n, numel(genes));
model = struct();
model.id = 'reporter';
model.rxns = cell(n, 1);
model.lb = zeros(n, 1);
model.ub = zeros(n, 1);
model.grRules = cell(n, 1);
for i = 1:n
    r = reactions{i};
    model.rxns{i} = char(r.id);
    model.lb(i) = double(r.lb);
    model.ub(i) = double(r.ub);
    model.grRules{i} = char(r.rule);
    for m = as_cellstr(r.substrates)
        S(strcmp(mets, m{1}), i) = -1;
    end
    for m = as_cellstr(r.products)
        S(strcmp(mets, m{1}), i) = 1;
    end
    for g = rule_genes(char(r.rule))
        rxnGeneMat(i, strcmp(genes, g{1})) = 1;
    end
end
model.mets = mets(:);
model.metNames = mets(:);
model.S = sparse(S);
model.rev = double(model.lb < 0);
model.c = zeros(n, 1);
model.b = zeros(numel(mets), 1);
model.comps = {'c'};
model.compNames = {'c'};
model.metComps = ones(numel(mets), 1);
model.rxnNames = model.rxns;
model.genes = genes(:);
model.rxnGeneMat = sparse(rxnGeneMat);
end


function out = rule_genes(rule)
tokens = regexp(rule, '[A-Za-z0-9_]+', 'match');
out = tokens(~ismember(tokens, {'and', 'or'}));
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
