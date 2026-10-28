function results = find_potential_errors_toy(ctx)
% MATLAB side of the findPotentialErrors scenario.

rules = ctx.inputs.rules;
if ischar(rules)
    rules = {rules};
end
rules = rules(:);
n = numel(rules);
model = struct();
model.rxns = arrayfun(@(i) sprintf('R%d', i), (1:n)', 'UniformOutput', false);
model.grRules = cellfun(@char, rules, 'UniformOutput', false);
issues = findPotentialErrors(model);
flagged = {};
for i = 1:numel(issues)
    flagged{end+1, 1} = issues(i).rxn; %#ok<AGROW>
end
flagged = sort(flagged);
if isempty(flagged)
    flagged = {};
end
results = struct('flagged', {flagged});
end
