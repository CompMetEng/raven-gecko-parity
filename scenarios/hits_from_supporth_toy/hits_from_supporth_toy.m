function results = hits_from_supporth_toy(ctx)
% MATLAB side of the hitsFromSuppOrth scenario.

inputs = ctx.inputs;
cases = as_list(inputs.cases);
results = struct();
for k = 1:numel(cases)
    c = cases{k};
    args = {'minSupport', double(c.min_support), 'requireTools', as_cellstr(c.require_tools)};
    if ~isempty(c.min_identity)
        args = [args, {'minIdentity', double(c.min_identity)}]; %#ok<AGROW>
    end
    if ~isempty(c.max_evalue)
        args = [args, {'maxEvalue', double(c.max_evalue)}]; %#ok<AGROW>
    end
    if logical(c.use_gene_map)
        args = [args, {'proteinToGene', char(inputs.gene_map)}]; %#ok<AGROW>
    end
    if ~isempty(c.species_to_model)
        names = fieldnames(c.species_to_model);
        pairs = cell(numel(names), 2);
        for i = 1:numel(names)
            pairs{i, 1} = names{i};
            value = c.species_to_model.(names{i});
            if isempty(value)
                value = '';
            end
            pairs{i, 2} = char(value);
        end
        args = [args, {'speciesToModel', pairs}]; %#ok<AGROW>
    end
    state = warning('off', 'all');
    try
        [bs, info] = hitsFromSuppOrth(char(inputs.table), char(c.organism), args{:});
    catch
        warning(state);
        results.(char(c.name)) = struct('error', true);
        continue
    end
    warning(state);
    pairs = {};
    for e = 1:numel(bs)
        for q = 1:numel(bs(e).fromGenes)
            pairs{end+1} = sprintf('%s|%s|%s|%s|%s|%g|%g', bs(e).fromId, bs(e).toId, bs(e).fromGenes{q}, bs(e).toGenes{q}, ...
                sprintf('%.10g', bs(e).bitscore(q)), bs(e).evalue(q), bs(e).identity(q)); %#ok<AGROW>
        end
    end
    support = cell(1, height(info));
    for q = 1:height(info)
        support{q} = sprintf('%s|%s|%s|%d|%s', info.modelId{q}, info.newGene{q}, info.templateGene{q}, info.support(q), info.supportedBy{q});
    end
    results.(char(c.name)) = struct('pairs', {sort(pairs)}, 'support', {sort(support)});
end
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
if iscell(value)
    out = value(:)';
else
    out = num2cell(value(:)');
end
end
