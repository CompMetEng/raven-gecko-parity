% Item 1: replace_metabolite_toy under replaceMets / contractModel
addpath(genpath([getenv('WORK') '/RAVEN']));
function model = nm()
model.mets = {'a';'b';'x';'d'};
model.metNames = {'oxygen';'o2';'x';'d'};
model.metComps = [1;1;1;1];
model.comps = {'c'};
model.compNames = {'cytoplasm'};
model.rxns = {'r1';'r2';'unrelated1';'unrelated2'};
model.rxnNames = model.rxns;
model.S = sparse([-1 0 0 0; 0 -1 0 0; 1 1 0 0; 0 0 -1 -1]);
model.lb = -1000*ones(4,1); model.ub = 1000*ones(4,1); model.rev = ones(4,1);
model.c = zeros(4,1); model.b = zeros(4,1);
end
function show(tag, m)
fprintf('--- %s\nmets: %s\nrxns: %s\nS=\n', tag, strjoin(m.mets','|'), strjoin(m.rxns','|')); disp(full(m.S));
end
m = nm(); show('(a) before replaceMets', m);
r = replaceMets(m,'oxygen','o2'); show('(b) after replaceMets(model,oxygen,o2)', r);
% (c) by hand
model = nm();
metabolite='oxygen'; replacement='o2';
repIdx = find(strcmp(replacement,model.metNames)); repMetIds = model.mets(repIdx);
metIdx = find(strcmp(metabolite,model.metNames));
rxnsWithMet = find(model.S(metIdx,:));
model.metNames(metIdx) = model.metNames(repIdx(1));
metCompsN = cellstr(num2str(model.metComps));
map = containers.Map(cellstr(num2str(transpose(1:length(model.comps)))),model.comps);
metCompsN = map.values(metCompsN);
metCompsN = strcat(lower(model.metNames),'[',metCompsN,']');
idxDelete=[];
for i = 1:length(repIdx)
    k=find(strcmp(metCompsN(repIdx(i)), metCompsN));
    if length(k)>1
        for j=2:length(k)
            model.S(k(1),:) = model.S(k(1),:)+model.S(k(j),:); idxDelete=[idxDelete;k(j)];
        end
    end
end
fprintf('--- (c) by hand: repMetIds = {%s}, idxDelete = %s\n', strjoin(repMetIds',','), mat2str(idxDelete));
model.S(idxDelete,:)=[]; model.mets(idxDelete)=[]; model.metNames(idxDelete)=[]; model.metComps(idxDelete)=[]; model.b(idxDelete,:)=[];
show('(c) after deletion step, before contractModel', model);
repIdx2 = find(ismember(model.mets, repMetIds));
fprintf('repIdx (as in replaceMets last lines) = %s  (isempty=%d)\n', mat2str(repIdx2), isempty(repIdx2));
modelS = sortModel(model);
x = [modelS.S; model.rev']';
disp('x before extra column:'); disp(full(x));
if ~isempty(repIdx2)
    toCheck = getIndexes(model,repIdx2,'mets'); [~, toCheck] = find(modelS.S(toCheck,:));
    x = [x, (1:size(x,1))']; x(toCheck,end)=0;
    fprintf('toCheck = %s\n', mat2str(toCheck')); disp('x with extra column:'); disp(full(x));
else
    disp('repIdx empty: contractModel(model,''mets'',[]) takes the no-restriction branch (p.mets empty)');
end
[~,I,J] = unique(x,'rows','first'); fprintf('unique rows: I=%s J=%s\n', mat2str(I'), mat2str(J'));
c1 = contractModel(model,'mets',repIdx2); show('contractModel(model,mets,repIdx)', c1);
c2 = contractModel(model); show('(d) contractModel(model) without mets', c2);
% what if repIdx pointed at the surviving metabolite (b was merged into a)
c3 = contractModel(model,'mets',find(strcmp(model.mets,'a'))); show('contractModel(model,mets,idx of surviving a)', c3);
