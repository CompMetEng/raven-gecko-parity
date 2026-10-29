classdef HumanGEMAdapterParity < TomlAdapter
% HumanGEMAdapterParity  The light_ecModel tutorial adapter, redirected to this scenario's PhylDist.mat.
%
% The parameters are those of the tutorial's model_adapter.toml, the file geckopy reads for
% the same project. TomlAdapter does not override getPhylDistStructPath, so it would use
% RAVEN's external/kegg/keggPhylDist.mat, a KEGG taxonomy artefact that is not vendored; see
% YeastGEMAdapterParity in ec_model_full_yeastgem for the reasoning behind the two-organism
% stand-in.
    methods
        function obj = HumanGEMAdapterParity()
            obj@TomlAdapter(fullfile(findGECKOroot, 'tutorials', 'light_ecModel'));
        end

        function x = getPhylDistStructPath(obj)
            [scenarioDir, ~, ~] = fileparts(mfilename('fullpath'));
            x = fullfile(scenarioDir, 'PhylDist.mat');
        end
    end
end
