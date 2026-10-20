% Item 4: probes P1-P9 on RAVEN develop
W = getenv('WORK'); addpath(genpath([W '/RAVEN']));
setRavenSolver('glpk');
sy = fullfile(W,'RAVEN','tutorials','data','models','smallYeast.yml');
tmp = tempname; mkdir(tmp);
runp(1,@() p1(sy)); runp(2,@() p2(sy)); runp(3,@() p3(sy)); runp(4,@() p4f());
runp(5,@() p5(W,tmp)); runp(6,@() p6(sy,tmp)); runp(7,@() p7(sy,tmp)); runp(8,@() p8()); runp(9,@() p9(tmp));
exit(0)

function runp(n,f)
try, f();
catch e, st=e.stack; loc=''; if ~isempty(st), loc=sprintf('%s:%d',st(1).name,st(1).line); end; fprintf('PROBE %d ERROR=%s | %s | at %s\n', n, e.identifier, strrep(strtrim(e.message),newline,' / '), loc); end
end
function pr(n,k,v)
if isnumeric(v)||islogical(v), v=mat2str(v); elseif iscell(v), v=['{' strjoin(cellfun(@char,v(:)','UniformOutput',false),' ; ') '}']; end
fprintf('PROBE %d %s=%s\n', n, k, v);
end

function p1(sy)
m1 = readYAMLmodel(sy); m2 = readYAMLmodel(sy);
nm = numel(m1.mets);
orph.mets={'ORPH_m'}; orph.metNames={'orphanMet'}; orph.compartments={m1.comps{1}};
m1 = addMets(m1, orph);
gs.genes={'ORPH_g'}; m1 = addGenesRaven(m1, gs);
ci = find(m1.ub>0,1); m1.c(:)=0; m1.c(ci)=1;
m2.id = 'second';
r.rxns={'NEWRXN'}; r.mets={m2.mets(1:2)'}; r.stoichCoeffs={[-1 1]}; r.lb=0; r.ub=1000;
m2 = addRxns(m2, r, 1);
pr(1,'model1_mets',numel(m1.mets)); pr(1,'model1_genes',numel(m1.genes)); pr(1,'model1_nnz_c',nnz(m1.c));
g = mergeModels({m1,m2});
pr(1,'merged_mets',numel(g.mets)); pr(1,'merged_genes',numel(g.genes)); pr(1,'merged_nnz_c',nnz(g.c));
pr(1,'orphan_met_survives',any(strcmp(g.mets,'ORPH_m'))); pr(1,'orphan_gene_survives',any(strcmp(g.genes,'ORPH_g')));
pr(1,'merged_rxns',numel(g.rxns)); pr(1,'NEWRXN_present',any(strcmp(g.rxns,'NEWRXN')));
end

function p2(sy)
m = readYAMLmodel(sy);
single = ~contains(m.grRules,{' and ',' or '}) & ~cellfun(@isempty,m.grRules);
cnt = full(sum(m.rxnGeneMat~=0,1)); cand = [];
for gi=1:numel(m.genes)
    if cnt(gi)==1 && any(single & m.rxnGeneMat(:,gi)~=0), cand=gi; break; end
end
g = m.genes{cand};
unusedM = @(x) sum(~any(x.S~=0,2)); unusedG = @(x) sum(~any(x.rxnGeneMat~=0,1));
pr(2,'gene',g); pr(2,'rxns_with_gene',m.rxns(m.rxnGeneMat(:,cand)~=0));
pr(2,'before_mets',numel(m.mets)); pr(2,'before_genes',numel(m.genes)); pr(2,'before_rxns',numel(m.rxns));
pr(2,'before_unused_mets',unusedM(m)); pr(2,'before_unused_genes',unusedG(m));
r = removeGenes(m, {g});
pr(2,'after_mets',numel(r.mets)); pr(2,'after_genes',numel(r.genes)); pr(2,'after_rxns',numel(r.rxns));
pr(2,'after_unused_mets',unusedM(r)); pr(2,'after_unused_genes',unusedG(r));
pr(2,'grRule_PGI_before',m.grRules(strcmp(m.rxns,'PGI'))); pr(2,'grRule_PGI_after',r.grRules(strcmp(r.rxns,'PGI')));
pr(2,'removed_gene_still_in_genes',any(strcmp(r.genes,g)));
pr(2,'rxn_still_present',any(ismember(m.rxns(m.rxnGeneMat(:,cand)~=0), r.rxns)));
end

function p3(sy)
m = readYAMLmodel(sy);
cidx = find(strcmp(m.comps,'c')); eidx = find(strcmp(m.comps,'m')); % smallYeast has c and m, no e
pr(3,'comps',m.comps);
nc = m.metNames(m.metComps==cidx); ne = m.metNames(m.metComps==eidx);
both = intersect(nc, ne);
pr(3,'n_names_in_both',numel(both)); 
% pick a name in both c and e that has no transport reaction yet
pick = both{1};
[m2, added] = addTransport(m,'c',{'m'},'metNames',{pick});
pr(3,'metName',pick); pr(3,'added_rxns',added);
[~,ia] = ismember(added, m2.rxns);
if isfield(m2,'subSystems')
    ss = m2.subSystems(ia);
    pr(3,'class_subSystems_entries',class(ss{1}));
    for k=1:numel(ss), v=ss{k}; if iscell(v), v=['{' strjoin(v(:)',' ; ') '}']; end; pr(3,sprintf('subSystems[%d]',k),v); end
else, pr(3,'subSystems_field','absent'); end
pr(3,'names_of_added',m2.rxnNames(ia));
end

function p4f()
model.genes = {'G1';'G2';'G3';'G4';'G5';'G6';'G7'};
model.grRules = {'G2 or G1';'G4 or G3 or G5';'G7 or G6'};
sc = [-1;-1;2;1;-3;1;1];
[nm, rem] = removeLowScoreGenes(model, sc);
pr(4,'grRules',nm.grRules); pr(4,'removed_genes',rem); pr(4,'genes_left',nm.genes);
end

function p5(W,tmp)
src = fileread(fullfile(W,'RAVEN','tests','function_tests','test_data','test_tasks.txt'));
lines = splitlines(src); hdr = lines{1};
cols = strsplit(hdr, char(9));
pr(5,'header_columns',cols);
ix = @(name) find(strcmp(cols,name));
mk = @(v) sprintf('%s', strjoin(v, char(9)));
base = repmat({''},1,numel(cols));
% (a) SHOULD FAIL variants
variants = {'0','"0"','1',''};
for k=1:numel(variants)
    row = mkrow(base, ix, 'ID','T1','DESCRIPTION','d','SHOULD FAIL',variants{k},'IN','A[c]','IN LB','0','IN UB','1000','OUT','B[c]','OUT LB','0','OUT UB','1000');
    f = fullfile(tmp, sprintf('a%d.txt',k)); fid=fopen(f,'w'); fprintf(fid,'%s\n%s\n',hdr,strjoin(row,char(9))); fclose(fid);
    try, ts = parseTaskList(f); pr(5,sprintf('a_shouldFail_input[%s]',variants{k}), ts(1).shouldFail);
    catch e, pr(5,sprintf('a_input[%s]_ERROR',variants{k}), strrep(e.message,newline,' / ')); end
end
% (b) IN separators
ins = {'A[c];;B[c]','A[c]; B[c]'};
for k=1:numel(ins)
    row = mkrow(base, ix, 'ID','T1','DESCRIPTION','d','IN',ins{k},'IN LB','0','IN UB','1000','OUT','C[c]','OUT LB','0','OUT UB','1000');
    f = fullfile(tmp, sprintf('b%d.txt',k)); fid=fopen(f,'w'); fprintf(fid,'%s\n%s\n',hdr,strjoin(row,char(9))); fclose(fid);
    try, ts = parseTaskList(f); pr(5,sprintf('b_IN[%s]_numel_inputs',ins{k}), numel(ts(1).inputs)); pr(5,sprintf('b_IN[%s]_inputs',ins{k}), ts(1).inputs);
    catch e, pr(5,sprintf('b_IN[%s]_ERROR',ins{k}), strrep(e.message,newline,' / ')); end
end
% (c) CHANGED RXN with empty LB/UB
row = mkrow(base, ix, 'ID','T1','DESCRIPTION','d','IN','A[c]','IN LB','0','IN UB','1000','OUT','C[c]','OUT LB','0','OUT UB','1000','CHANGED RXN','R1');
f = fullfile(tmp,'c.txt'); fid=fopen(f,'w'); fprintf(fid,'%s\n%s\n',hdr,strjoin(row,char(9))); fclose(fid);
try, ts = parseTaskList(f); pr(5,'c_changed',ts(1).changed); pr(5,'c_LBrxn',ts(1).LBrxn); pr(5,'c_UBrxn',ts(1).UBrxn);
catch e, pr(5,'c_ERROR',strrep(e.message,newline,' / ')); end
% (c') the same on a continuation row of a task
row0 = mkrow(base, ix, 'ID','T1','DESCRIPTION','d','IN','A[c]','IN LB','0','IN UB','1000','OUT','C[c]','OUT LB','0','OUT UB','1000');
row1 = mkrow(base, ix, 'CHANGED RXN','R1');
f = fullfile(tmp,'c2.txt'); fid=fopen(f,'w'); fprintf(fid,'%s\n%s\n%s\n',hdr,strjoin(row0,char(9)),strjoin(row1,char(9))); fclose(fid);
try, ts = parseTaskList(f); pr(5,'c2_changed',ts(1).changed); pr(5,'c2_LBrxn',ts(1).LBrxn); pr(5,'c2_UBrxn',ts(1).UBrxn);
catch e, pr(5,'c2_ERROR',strrep(e.message,newline,' / ')); end
end

function row = mkrow(base, ix, varargin)
row = base; for k=1:2:numel(varargin), row{ix(varargin{k})} = varargin{k+1}; end
end

function p6(sy,tmp)
m = readYAMLmodel(sy);
m0 = m; m0.metNames{3} = '';
try, writeYAMLmodel(m0, fullfile(tmp,'probe6_empty.yml')); pr(6,'empty_metName_write','succeeded');
catch e, pr(6,'empty_metName_write_ERROR',[e.identifier ' | ' strrep(strtrim(e.message),newline,' / ')]); end
m.lb(2) = 1/3; if m.ub(2) < m.lb(2), m.ub(2)=1; end
f = fullfile(tmp,'probe6.yml'); writeYAMLmodel(m, f);
L = splitlines(fileread(f));
rid = m.rxns{2}; kr = find(contains(L, ['id: ' rid]),1);
pr(6,'reaction_block',strjoin(L(kr:min(kr+8,numel(L)))',' / '));
kc = find(contains(L,'compartments'),1);
if isempty(kc), pr(6,'compartments_block','no line contains compartments'); else, pr(6,'compartments_block',strjoin(L(kc:min(kc+2*numel(m.comps)+1,numel(L)))',' / ')); end
pr(6,'model_comps_order',m.comps);
end

function p7(sy,tmp)
m = readYAMLmodel(sy);
d = fullfile(tmp,'git7'); mkdir(d);
cwd = pwd; cleanup = onCleanup(@() cd(cwd)); cd(d);
exportForGit(m,'prefix','probe','path',d);
fl = dir(fullfile(d,'**','*')); fl = fl(~[fl.isdir]);
pr(7,'files',arrayfun(@(x) fullfile(extractAfter(x.folder,d), x.name), fl,'UniformOutput',false));
[~,~,ext] = cellfun(@fileparts,{fl.name},'UniformOutput',false);
pr(7,'extensions',unique(ext));
end

function p8()
rules = {'(A and B) or C','A && B','(A or B) and C','a OR b','A and',''};
for k=1:numel(rules)
    try, pr(8,sprintf('isDnfGrRule[%s]',rules{k}), isDnfGrRule(rules{k}));
    catch e, pr(8,sprintf('isDnfGrRule[%s]_ERROR',rules{k}), strrep(e.message,newline,' / ')); end
end
end

function p9(tmp)
m.id='dup'; m.mets={'a1';'a2';'b'}; m.metNames={'A';'A';'B'}; m.metComps=[1;1;1]; m.comps={'c'}; m.compNames={'cytosol'};
m.rxns={'a2_to_b'}; m.rxnNames=m.rxns; m.S=sparse([0;-1;1]); m.lb=0; m.ub=1000; m.rev=0; m.c=0; m.b=zeros(3,1);
m.genes={}; m.grRules={''}; m.rxnGeneMat=sparse(1,0);
hdr = strjoin({'','ID','DESCRIPTION','IN','IN LB','IN UB','OUT','OUT LB','OUT UB'},char(9));
row = strjoin({'','T1','dup names','A[c]','1','1000','B[c]','1','1000'},char(9));
f = fullfile(tmp,'p9.txt'); fid=fopen(f,'w'); fprintf(fid,'%s\n%s\n',hdr,row); fclose(fid);
[rep, ess, ts] = checkTasks(m, f, false, false, true);
pr(9,'taskReport_ok',rep.ok); pr(9,'taskReport_description',rep.description); pr(9,'essentialRxns',ess);
end
