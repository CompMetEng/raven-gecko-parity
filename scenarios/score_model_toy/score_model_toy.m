function results = score_model_toy(ctx)
% MATLAB side of the scoreModel scenario with zero and missing expression levels.

inputs = ctx.inputs;
model = readYAMLmodel(inputs.model);
genes = sort(model.genes(:));
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    arrayData = struct();
    arrayData.genes = genes(1:6);
    arrayData.tissues = {'t'};
    arrayData.levels = double(c.levels(:));
    arrayData.threshold = double(c.threshold);
    [rxnScores, geneScores] = scoreModel(model, [], arrayData, 't', ...
        'noGeneScore', double(inputs.no_gene_score), ...
        'isozymeScoring', char(inputs.isozyme_scoring), ...
        'complexScoring', char(inputs.complex_scoring));

    records = cell(1, 6);
    for g = 1:6
        i = find(strcmp(model.genes, genes{g}), 1);
        hasScore = ~isnan(geneScores(i));
        score = 0;
        if hasScore
            score = double(geneScores(i));
        end
        records{g} = struct('gene', genes{g}, 'has_score', hasScore, 'score', score);
    end
    [sortedRxns, order] = sort(model.rxns);
    reactions = cell(1, numel(sortedRxns));
    for r = 1:numel(sortedRxns)
        reactions{r} = struct('reaction', sortedRxns{r}, 'score', double(rxnScores(order(r))));
    end
    results.(char(c.name)) = struct('gene_scores', {records}, 'reaction_scores', {reactions});
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
