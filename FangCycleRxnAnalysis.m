load('MTBmodels.mat')
load('BJmap.mat')
load('Bigg2Kegg.mat')
initCobraToolbox

%% Find the unconstrained reactions in Fang
[v1F v2F] = fluxVariability(iNJ661m);
xF = find(all(abs([v1F v2F]) > 999.99,2));
x1F = find(v1F < -999.99);
x2F = find(v2F > 999.99);
x1Finfo = [v1F(x1F) v2F(x1F) iNJ661m.lb(x1F) iNJ661m.ub(x1F)];
x2Finfo = [v1F(x2F) v2F(x2F) iNJ661m.lb(x2F) iNJ661m.ub(x2F)];
xFinfo = printRxnFormula(iNJ661m,iNJ661m.rxns(xF));
subsysxF = iNJ661m.subSystems(xF);
genF = findGenesFromRxns(iNJ661m,iNJ661m.rxns(xF));
xFinfo = [iNJ661m.rxns(xF) xFinfo genF subsysxF];
rxnInfo(iNJ661m,xF)

%% Find the unconstrained reactions in Jamshidi
[v1J v2J] = fluxVariability(jamshidiBiGG);
xJ = find(all(abs([v1J v2J]) > 999.99,2));
x1J = find(v1J < -999.99);
x2J = find(v2J > 999.99);
x1Jinfo = [v1J(x1J) v2J(x1J) jamshidiBiGG.lb(x1J) jamshidiBiGG.ub(x1J)];
x2Jinfo = [v1J(x2J) v2J(x2J) jamshidiBiGG.lb(x2J) jamshidiBiGG.ub(x2J)];
xJinfo = printRxnFormula(jamshidiBiGG,jamshidiBiGG.rxns(xJ));
xJinfo = [jamshidiBiGG.rxns(xJ) xJinfo];
rxnInfo(jamshidiBiGG,xJ)

%% Common reactions in Fang and Jamshidi
[c iF iJ] = intersect(iNJ661m.rxns,jamshidiBiGG.rxns);

%% The neighboring reactions of the unconstrained reactions
neighborinfoF = cell(size(xF,1),3);
for i = 1:size(xF,1)
    [nr1 ng1 m] = findNeighborRxns(iNJ661m,iNJ661m.rxns(xF(i)));
    neighborinfoF{i,1} = nr1;
    neighborinfoF{i,2} = ng1;
    neighborinfoF{i,3} = m;
end
Fangneighborinfo = cell(20,1);
for i = 1:20
    Fangneighborinfo{i,1} = [neighborinfoF{i,3} neighborinfoF{i,1}' neighborinfoF{i,2}'];
end

%% Extracting the subsystems which have the unconstrained reactions
PurineFang = extractSubSysModel(iNJ661m,subsysxF{1});
NucSugarFang = extractSubSysModel(iNJ661m,subsysxF{2});
TCAFang = extractSubSysModel(iNJ661m,subsysxF{4});
RedoxFang = extractSubSysModel(iNJ661m,subsysxF{5});
GluFang = extractSubSysModel(iNJ661m,subsysxF{7});
TransportFang = extractSubSysModel(iNJ661m,subsysxF{8});
PyrimidineFang = extractSubSysModel(iNJ661m,subsysxF{14});
UbiquinoneFang = extractSubSysModel(iNJ661m,subsysxF{16});
CoAFang = extractSubSysModel(iNJ661m,subsysxF{20});

CbModeltoEXPA(PurineFang,'PurineFang.expa')
CbModeltoEXPA(NucSugar,'NucSugarFang.expa')
CbModeltoEXPA(NucSugarFang,'NucSugarFang.expa')
CbModeltoEXPA(TCAFang,'TCAFang.expa')
CbModeltoEXPA(RedoxFang,'RedoxFang.expa')
CbModeltoEXPA(GluFang,'GluFang.expa')
CbModeltoEXPA(TransportFang,'TransportFang.expa')
CbModeltoEXPA(PyrimidineFang,'PyrimidineFang.expa')
CbModeltoEXPA(UbiquinoneFang,'UbiquinoneFang.expa')
CbModeltoEXPA(CoAFang,'CoAFang.expa')

%% Find the reactions mapped to each key metabolite

[FangATPrxns] = findRxnsFromMets(iNJ661m,'atp[c]');
[FangNADrxns] = findRxnsFromMets(iNJ661m,'nad[c]');
[FangNADHrxns] = findRxnsFromMets(iNJ661m,'nadh[c]');
[FangNADPHrxns] = findRxnsFromMets(iNJ661m,'naph[c]');
[FangNADPrxns] = findRxnsFromMets(iNJ661m,'nap[c]');
[FangADPrxns] = findRxnsFromMets(iNJ661m,'adp[c]');

[FangATPerxns] = findRxnsFromMets(iNJ661m,'atp[e]');

ex = strmatch('EX_',iNJ661m.rxns);

% Extract subnetwork of all reactions involving ATP
ATPFang = extractSubNetwork(iNJ661m,[FangATPrxns;FangADPrxns;FangATPerxns;iNJ661m.rxns(ex)]);
CbModeltoEXPA(ATPFang,'ATPFang.expa')

%% Creating a GIMME model of Fang based on aerobic culture condition data
load('mtbPROMmodelInputs.mat', 'SG0')
load('mtbPROMmodelInputs.mat', 'SG1', 'SG2', 'SG3', 'SG5', 'SG7')
load('mtbPROMmodelInputs.mat', 'SGgen')
load('AlldataRv.mat', 'AlldataRv', 'Allnegcont')
SG0 = cell2mat(SG0);
for i = [1 2 3 5 7]
    eval(['SG' num2str(i) '=cell2mat(SG' num2str(i) ');']);
end
save mtbPROMmodelInputs SG* -append
SG0exp.Data = mean(SG0,2) > max(Allnegcont);
[x x1 x2] = intersect(iNJ661m.genes,SGgen);
SG0exp.Locus = x1;
SG0exp.Data = SG0exp.Data(x2);

[ix ix ix] = cellfun(@(x) intersect(x,SGgen),genF,'UniformOutput',false);
ix(:,2) = cellfun(@(x) SG0(x,:),ix(:,1),'UniformOutput',false);
ix(:,3) = cellfun(@(x) mean(x,2),ix(:,2),'UniformOutput',false);

medianRv = median(cell2mat(AlldataRv(:,2:end)),2);
ix(:,4) = cellfun(@(x) medianRv(x),ix(:,1),'UniformOutput',false);

medianControl = median(Allnegcont)
control95 = quantile(Allnegcont,0.95)
control33 = quantile(Allnegcont,0.33)
control99 = quantile(Allnegcont,0.99)

% [v1F v2F] = fluxVariability(iNJ661m);
xFbig = find(all(abs([v1F v2F]) > 999.99,2));
xFsmall = find(all(abs([v1F v2F]) < 10^-8,2));
lowfluxrxnsF = [iNJ661m.rxns(xFsmall) iNJ661m.subSystems(xFsmall)];

save MTBConditionSpecifModel iNJ661m SG0exp
[SG0FangGIMME SG0FangGIMMErxns] = createTissueSpecificModelSM(iNJ661m,SG0exp);
% [SG0FangIMAT SG0FangIMATrxns] = createTissueSpecificModelSM(iNJ661m,SG0exp,1,1,[],'iMAT');

[v1FG v2FG] = fluxVariability(SG0FangGIMME);
xFG = find(all(abs([v1FG v2FG]) > 999.99,2));
SG0FangGIMME.rxns(xFG)
xFGbig = find(all(abs([v1FG v2FG]) > 999.99,2));
xFGsmall = find(all(abs([v1FG v2FG]) < 10^-8,2));
lowfluxrxnsFG = [SG0FangGIMME.rxns(xFGsmall) SG0FangGIMME.subSystems(xFGsmall)];

[FangGIMMEred GFhasflux FGvmax FGvmin] = reduceModel(SG0FangGIMME,10^-10);
soln = optimizeCbModel(FangGIMMEred)
soln = optimizeCbModel(SG0FangGIMME)

CbModeltoEXPA(FangGIMMEred,'FangGIMMEreduced.expa')
[FangGIMMEred8 GF8hasflux FGvmax FGvmin] = reduceModel(SG0FangGIMME,10^-8);
soln = optimizeCbModel(FangGIMMEred8)
isSameCobraModel(FangGIMMEred8,FangGIMMEred)
[FangGIMMEred7 GF7hasflux FGvmax FGvmin] = reduceModel(SG0FangGIMME,10^-7);
soln = optimizeCbModel(FangGIMMEred7)
[FangGIMMEred6 GF6hasflux FGvmax FGvmin] = reduceModel(SG0FangGIMME,10^-6);
CbModeltoEXPA(FangGIMMEred6,'FangGIMMEreduced6.expa')
save FangRxnInfo xF* *F *G SG0FangGIMME*
save FangRxnInfo Fang* -append
%% 

%% Elementary Flux Mode stuff
load('test_out.mat', 'E')

load('FangEM1000.mat')
load('iNJ661mCorrRxnSets.mat', 'setsSorted')
load('iNJ661mCorrRxnSets.mat', 'setNoSorted')
load('MTBmodels.mat', 'iNJ661m')

EFMrxns = cell(31,2);for i=1:size(E,1);EFMrxns{i,1} = find(E(i,:) ~= 0);EFMrxns{i,2} = E(i,EFMrxns{i,1});end
EFMrxns(:,3) = cellfun(@(x) iNJ661m.rxns(x),EFMrxns(:,1),'UniformOutput',false);

xFbigSet = setNoSorted(xF);
xFbigSet = num2cell(xFbigSet);
for i = 1:20; if xFbigSet{i,1} ~= 0; xFbigSet{i,2} = setsSorted{xFbigSet{i,1}}.names;end;end
load('FangRxnInfo.mat', 'Fangneighborinfo')
save FangRxnInfo xFbigSet -append