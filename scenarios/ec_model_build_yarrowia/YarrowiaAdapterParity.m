classdef YarrowiaAdapterParity < TomlAdapter
% YarrowiaAdapterParity  The scenario's model_adapter.toml, with the eciYali5-GEM data folder and this scenario's PhylDist.mat.
%
% The adapter folder of the scenario holds only the TOML file; the data files and the
% conventional model live in a checkout of eciYali5-GEM whose folder is given to the
% constructor, and the same two paths are set on the geckopy adapter. TomlAdapter does not
% override getPhylDistStructPath, so the two-organism stand-in of ec_model_full_yeastgem is
% used instead of RAVEN's KEGG taxonomy file, which is not vendored.
    methods
        function obj = YarrowiaAdapterParity(eciYaliDir)
            [scenarioDir, ~, ~] = fileparts(mfilename('fullpath'));
            obj@TomlAdapter(scenarioDir);
            obj.params.path = char(eciYaliDir);
            obj.params.convGEM = fullfile(char(eciYaliDir), 'model', 'iYali5-GEM', 'iYali5.xml');
        end

        function x = getPhylDistStructPath(obj)
            [scenarioDir, ~, ~] = fileparts(mfilename('fullpath'));
            x = fullfile(scenarioDir, 'PhylDist.mat');
        end
    end
end
