function results = hpa_scores_toy(ctx)
% MATLAB side of the HPA branch of the scoreModel scenario.

inputs = ctx.inputs;
model = readYAMLmodel(inputs.model);
genes = sort(model.genes(:));
genes = genes(1:6);
baseRows = as_list(inputs.rows);
cases = as_list(inputs.cases);

results = struct();
for k = 1:numel(cases)
    c = cases{k};
    rows = baseRows;
    if isfield(c, 'extra_rows')
        rows = [rows, as_list(c.extra_rows)];
    end
    hpaData = build_hpa(rows, genes);
    multiple = 'average';
    if strcmp(char(c.multiple), 'best')
        multiple = 'max';
    end
    args = {'multipleCellScoring', multiple};
    if ~isempty(char(c.celltype))
        args = [args, {'celltype', char(c.celltype)}]; %#ok<AGROW>
    end
    if isfield(c, 'table')
        entries = as_list(c.table);
        levelScores = struct('names', {cellfun(@(e) char(e.name), entries, 'UniformOutput', false)}, ...
            'scores', cellfun(@(e) double(e.score), entries));
        args = [args, {'hpaLevelScores', levelScores}]; %#ok<AGROW>
    end
    try
        [~, ~, hpaScores] = scoreModel(model, hpaData, [], char(c.tissue), args{:});
    catch
        results.(char(c.name)) = struct('error', true);
        continue;
    end
    records = cell(1, numel(genes));
    for g = 1:numel(genes)
        i = find(strcmp(model.genes, genes{g}), 1);
        has = isfinite(hpaScores(i));
        score = 0;
        if has
            score = double(hpaScores(i));
        end
        records{g} = struct('gene', genes{g}, 'has_score', has, 'score', score);
    end
    results.(char(c.name)) = records;
end
end


function hpaData = build_hpa(rows, genes)
nRows = numel(rows);
tissueOf = cell(nRows, 1); celltypeOf = cell(nRows, 1); levelOf = cell(nRows, 1); geneOf = zeros(nRows, 1);
for i = 1:nRows
    tissueOf{i} = char(rows{i}.tissue);
    celltypeOf{i} = char(rows{i}.celltype);
    levelOf{i} = char(rows{i}.level);
    geneOf(i) = double(rows{i}.gene) + 1;
end
[~, firstIdx] = unique(strcat(tissueOf, '|', celltypeOf), 'stable');
colKeys = strcat(tissueOf(firstIdx), '|', celltypeOf(firstIdx));
levels = unique(levelOf, 'stable');
gene2Level = zeros(numel(genes), numel(colKeys));
for i = 1:nRows
    col = find(strcmp(colKeys, [tissueOf{i} '|' celltypeOf{i}]));
    gene2Level(geneOf(i), col) = find(strcmp(levels, levelOf{i}));
end
hpaData = struct();
hpaData.genes = genes;
hpaData.tissues = tissueOf(firstIdx)';
hpaData.celltypes = celltypeOf(firstIdx)';
hpaData.levels = levels;
hpaData.gene2Level = sparse(gene2Level);
end


function out = as_list(value)
% jsondecode gives a struct array when every element has the same fields and a cell otherwise.
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end
