function varargout = pr_common(action, varargin)
% Helpers shared by the probes: JSON output with NaN as "NaN", error capture.
switch action
    case 'save'
        id = varargin{1}; s = pr_clean(varargin{2});
        fid = fopen(fullfile(getenv('PROBE_OUT'), [id '.json']), 'w');
        fwrite(fid, jsonencode(s, 'PrettyPrint', true), 'char'); fclose(fid);
    case 'err'
        e = varargin{1};
        varargout{1} = struct('identifier', e.identifier, 'message', e.message);
    case 'model'
        varargout{1} = readYAMLmodel(fullfile(getenv('RAVEN_ROOT'), 'tutorials', 'data', 'models', 'smallYeast.yml'));
    case 'unusedGenes'
        m = varargin{1}; out = {};
        for i = 1:numel(m.genes)
            pat = ['(^|[^A-Za-z0-9_.])' regexptranslate('escape', m.genes{i}) '($|[^A-Za-z0-9_.])'];
            if ~any(~cellfun(@isempty, regexp(m.grRules, pat, 'once')))
                out{end+1, 1} = m.genes{i}; %#ok<AGROW>
            end
        end
        varargout{1} = out;
    case 'stats'
        m = varargin{1};
        s = struct('genes', numel(m.genes), 'mets', numel(m.mets), 'rxns', numel(m.rxns));
        s.genesInNoRule = pr_common('unusedGenes', m);
        varargout{1} = s;
end
end

function x = pr_clean(x)
if isnumeric(x) && ~isempty(x) && any(isnan(x(:)))
    c = num2cell(x);
    c(isnan(x)) = {'NaN'};
    x = reshape(c, size(x));
elseif isstruct(x)
    for k = 1:numel(x)
        f = fieldnames(x);
        for j = 1:numel(f)
            x(k).(f{j}) = pr_clean(x(k).(f{j}));
        end
    end
elseif iscell(x)
    for k = 1:numel(x)
        x{k} = pr_clean(x{k});
    end
end
end
