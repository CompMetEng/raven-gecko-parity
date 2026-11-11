function results = okp_kcat_ectestgem(ctx)
% MATLAB side of the OpenKineticsPredictor (OKP) scenario, mirroring run.py.
%
% findMetSmiles, the job input that kcatsFromOKP writes (with the submission stopped before
% it reaches the network) and the stored result of a job read through kcatsFromOKP, on one
% ecModel built from ecTestGEM. See scenario.yml for scope (CSV boundary only, not the real OKP REST API) and for
% [GECKO #437], the real indexing bug this scenario's write checkpoint found,
% confirmed and fixed --- this scenario's own MATLAB result was generated against a
% worktree with that fix applied ahead of merge.
%
% findMetSmiles and kcatsFromOKP resolve fixture/output paths through the adapter's own
% params.path, unlike find_met_smiles, which takes an explicit path. A *copy* of the adapter
% --- ModelAdapter is a value class, so this cannot affect the original, same pattern as
% kcat_chain_ectestgem.m --- has its params.path redirected to this scenario's own data/
% folder for findMetSmiles, and to a scratch project for the two OKP checkpoints.

adapter = ModelAdapterManager.getAdapter(ctx.inputs.adapter_matlab);
model = loadConventionalGEM('modelAdapter', adapter);
ecModel = makeEcModel(model, false, adapter);
ecModel = getECfromGEM(ecModel);

okpAdapter = adapter;
okpAdapter.params.path = ctx.inputs.scenario_root;

[results.smiles, ecModel] = checkpoint_smiles(ecModel, okpAdapter);
% kcatsFromOKP resolves the API key after it has written the job input; without a key it
% stops there, so nothing is sent.
previousKey = getenv('OKP_API_KEY');
restoreKey = onCleanup(@() setenv('OKP_API_KEY', previousKey));
setenv('OKP_API_KEY', '');
[results.write, ecModel] = checkpoint_write(ecModel, adapter);
results.read_output = checkpoint_read(ecModel, adapter, ctx.inputs);

end


% --------------------------------------------------------------------------- %
% smiles: findMetSmiles
% --------------------------------------------------------------------------- %

function [out, ecModel] = checkpoint_smiles(ecModel, okpAdapter)
ecModel = findMetSmiles(ecModel, 'modelAdapter', okpAdapter, 'verbose', false);
out = smiles_snapshot(ecModel);
end


function out = smiles_snapshot(ecModel)
n = numel(ecModel.mets);
records = cell(1, n);
for k = 1:n
    s = '';
    if ~isempty(ecModel.metSmiles{k})
        s = ecModel.metSmiles{k};
    end
    records{k} = struct('met', ecModel.mets{k}, 'smiles', s);
end
[~, order] = sort(ecModel.mets);
records = records(order);
out = struct();
for k = 1:n
    out.(records{k}.met) = records{k}.smiles;
end
end


% --------------------------------------------------------------------------- %
% write: the job input of kcatsFromOKP
% --------------------------------------------------------------------------- %

function [out, ecModel] = checkpoint_write(ecModel, adapter)
% ec.rxns order on this fixture is [R2_EXP_1, R2_EXP_2, R2_REV_EXP_1, R2_REV_EXP_2,
% R3, R5] (confirmed by direct execution, identical on both sides). Requesting
% {R2_EXP_1, R2_EXP_2, R3} skips the two REV entries sitting in between --- exactly
% the subset shape [GECKO #437]'s bug mishandled, so this checkpoint doubles as its
% regression check.
withSmilesMask = ismember(ecModel.ec.rxns, {'R2_EXP_1', 'R2_EXP_2', 'R3'});
out.with_smiles = job_input(ecModel, adapter, withSmilesMask);

% m2 (R5's substrate) deliberately loses its SMILES here: kcatsFromOKP always drops an
% entry without a SMILES from the job input, so R5's entry does not appear at all.
ecModel.metSmiles(strcmp(ecModel.mets, 'm2c')) = {''};
fullMask = ismember(ecModel.ec.rxns, {'R2_EXP_1', 'R2_EXP_2', 'R3', 'R5'});
out.after_clearing_m2_with_filter = job_input(ecModel, adapter, fullMask);
end


function out = job_input(ecModel, adapter, mask)
% The pairs of the job input that kcatsFromOKP writes for the selected reactions, in a
% scratch project. The call fails at the API key, after the input is written.
project = tempname;
mkdir(fullfile(project, 'data'));
cleanup = onCleanup(@() rmdir(project, 's'));
scratch = adapter;
scratch.params.path = project;
try
    kcatsFromOKP(ecModel, 'method', 'CataPro', 'ecRxns', mask, 'save', false, ...
        'modelAdapter', scratch);
catch
end
csvFile = fullfile(project, 'data', 'OKP.csv');
fID = fopen(csvFile);
raw = textscan(fID, '%s %s', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fID);
pairs = strcat(raw{1}, ',', raw{2});
out = reshape(sort(pairs), 1, []);
end


% --------------------------------------------------------------------------- %
% read_output: the stored result of a job, read through kcatsFromOKP
% --------------------------------------------------------------------------- %

function out = checkpoint_read(ecModel, adapter, inputs)
project = tempname;
mkdir(fullfile(project, 'data'));
cleanup = onCleanup(@() rmdir(project, 's'));
copyfile(inputs.okp_result_path, fullfile(project, 'data', 'OKP_output_catapro.csv'));
fid = fopen(fullfile(project, 'data', 'OKP_job_catapro.txt'), 'w');
fprintf(fid, 'jobId: scenario-job\nmethod: CataPro\nresultOf: scenario-job\n');
fclose(fid);
scratch = adapter;
scratch.params.path = project;
[kcatList, done] = kcatsFromOKP(ecModel, 'method', 'CataPro', 'useStored', true, ...
    'save', false, 'modelAdapter', scratch);
n = numel(kcatList.rxns);
records = cell(1, n);
for k = 1:n
    % Recorded verbatim: 'OKP-CataPro', 'OKP-BRENDA' on both sides.
    records{k} = struct('reaction', kcatList.rxns{k}, 'gene', kcatList.genes{k}, ...
        'substrate', kcatList.substrates{k}, 'kcat', double(kcatList.kcats(k)), ...
        'source', kcatList.kcatSource{k});
end
keys = cellfun(@(r) [r.reaction char(1) r.gene], records, 'UniformOutput', false);
[~, order] = sort(keys);
out.done = logical(done);
out.rows = records(order);
end
