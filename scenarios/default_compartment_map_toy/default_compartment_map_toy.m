function results = default_compartment_map_toy(ctx) %#ok<INUSD>
% MATLAB side of the defaultCompartmentMap scenario.

map = defaultCompartmentMap();
labels = keys(map);
pairs = cell(1, numel(labels));
for i = 1:numel(labels)
    pairs{i} = sprintf('%s|%s', labels{i}, map(labels{i}));
end
results = struct('count', double(map.Count), 'pairs', {sort(pairs)});
end
