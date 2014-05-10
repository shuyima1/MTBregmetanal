initCobraToolbox
%% Metabolite-specific reaction analysis


load('MTBmodels.mat', 'iNJ661m')
ex = strmatch('EX_',iNJ661m.rxns);
Fangnoex = changeRxnBounds(iNJ661m,iNJ661m.rxns(ex),0,'b');
[v1nx v2nx] = fluxVariability(Fangnoex);

ix1 = find(abs(v1nx) > 0.01);
ix1(:,2) = v1nx(ix1);
ix2 = find(abs(v2nx) > 0.01);
ix2(:,2) = v2nx(ix2);

negrxns = iNJ661m.rxns(ix1(:,1));
posrxns = iNJ661m.rxns(ix2(:,1));
negrxns(:,2) = printRxnFormula(iNJ661m,negrxns(:,1));
posrxns(:,2) = printRxnFormula(iNJ661m,posrxns(:,1))

[FangnoexMets HiEnergyDemandRxns] = addDemandReaction(Fangnoex,iNJ661m.mets([133;181;312;313;390;402;445;535;536;537;538]))
FangnoexMetsATP = changeObjective(FangnoexMets,'DM_atp[c]');
solnATP = optimizeCbModel(FangnoexMetsATP)
FangnoexMetsNAD = changeObjective(FangnoexMets,'DM_nad[c]');
solnNAD = optimizeCbModel(FangnoexMetsNAD)
FangnoexMetsNADP = changeObjective(FangnoexMets,'DM_nadp[c]');
solnNADP = optimizeCbModel(FangnoexMetsNADP)
FangnoexMetsNADPH = changeObjective(FangnoexMets,'DM_nadph[c]');
solnNADPH = optimizeCbModel(FangnoexMetsNADPH)
FangnoexMetsNADH = changeObjective(FangnoexMets,'DM_nadh[c]');
solnNADH = optimizeCbModel(FangnoexMetsNADH)
FangnoexMetsADP = changeObjective(FangnoexMets,'DM_adp[c]');
solnADP = optimizeCbModel(FangnoexMetsADP)
FangnoexMetsFAD = changeObjective(FangnoexMets,'DM_fad[c]');
solnFAD = optimizeCbModel(FangnoexMetsFAD)
FangnoexMetsH = changeObjective(FangnoexMets,'DM_h[c]');
solnH = optimizeCbModel(FangnoexMetsH)
FangnoexMetsFADH2 = changeObjective(FangnoexMets,'DM_fadh2[c]');
solnFADH2 = optimizeCbModel(FangnoexMetsFADH2)
FangnoexMetsGTP = changeObjective(FangnoexMets,'DM_gtp[c]');
solnGTP = optimizeCbModel(FangnoexMetsGTP)
FangnoexMetsITP = changeObjective(FangnoexMets,'DM_itp[c]');
solnITP = optimizeCbModel(FangnoexMetsITP)

[v1nxATP v2nxATP] = fluxVariability(FangnoexMetsATP);[v1nxNAD v2nxNAD] = fluxVariability(FangnoexMetsNAD);[v1nxNADP v2nxNADP] = fluxVariability(FangnoexMetsNADP);[v1nxNADH v2nxNADH] = fluxVariability(FangnoexMetsNADH);[v1nxNADPH v2nxNADPH] = fluxVariability(FangnoexMetsNADPH);

ix1ATP = find(abs(v1nxATP) > 0.01);
ix1ATP(:,2) = v1nxATP(ix1ATP);
ix2ATP = find(abs(v2nxATP) > 0.01);
ix2ATP(:,2) = v2nxATP(ix2ATP);

negrxnsATP = iNJ661m.rxns(ix1ATP(:,1));
posrxnsATP = iNJ661m.rxns(ix2ATP(:,1));


cyclerxns = union(ix1(:,1),ix2(:,1));
cyclemets = [];
for i = 1:numel(cyclerxns)
    cyclemets=union(cyclemets,find(iNJ661m.S(:,cyclerxns(i)) ~= 0));
end

[FangnoexMets2 CMetDemandRxns]= addDemandReaction(Fangnoex,iNJ661m.mets(cyclemets));

x = strmatch('DM_',FangnoexMets2.rxns);
nzmet = [];
for i = 1:size(CMetDemandRxns,2)
    model = changeObjective(FangnoexMets2,CMetDemandRxns{i});
    soln = optimizeCbModel(model);
    if abs(soln.f) > 0.0001
        nzmet = [nzmet;i soln.f];
    end
end

atprxns = findRxnsFromMets(iNJ661m,'atp[c]');
atprxns(:,2) = cellfun(@(x) strmatch(x,iNJ661m.rxns,'exact'),atprxns(:,1),'UniformOutput',false);
atprxns(:,3) = cellfun(@(x) full(iNJ661m.S(181,x)),atprxns(:,2),'UniformOutput',false);

atpnocyclerxns = setdiff(atprxns,[posrxns(:,1);negrxns(:,1)]);
atpnocyclerxns(:,2) = cellfun(@(x) strmatch(x,iNJ661m.rxns,'exact'),atpnocyclerxns,'UniformOutput',false);
atpnocyclerxns(:,3) = cellfun(@(x) full(iNJ661m.S(181,x)),atpnocyclerxns(:,2),'UniformOutput',false);

atpflux = [v1nx(cell2mat(atpnocyclerxns(:,2))) v2nx(cell2mat(atpnocyclerxns(:,2)))];
sum(cell2mat(atpnocyclerxns(:,3)).*atpflux(:,1))
sum(cell2mat(atpnocyclerxns(:,3)).*atpflux(:,2))

atpflux2 = [v1nx(cell2mat(atprxns(:,2))) v2nx(cell2mat(atprxns(:,2)))];sum(cell2mat(atprxns(:,3)).*atpflux2(:,1))
sum(cell2mat(atprxns(:,3)).*atpflux2(:,2))

[v1 v2] = fluxVariability(iNJ661m);
atpfluxF = [v1(cell2mat(atpnocyclerxns(:,2))) v2(cell2mat(atpnocyclerxns(:,2)))];
sum(cell2mat(atpnocyclerxns(:,3)).*atpfluxF(:,1))
sum(cell2mat(atpnocyclerxns(:,3)).*atpfluxF(:,2))
sum(cell2mat(atpnocyclerxns(1:end-1,3)).*atpfluxF(1:end-1,1))
sum(cell2mat(atpnocyclerxns(1:end-1,3)).*atpfluxF(1:end-1,2))

load('PROMrun1.mat', 'vb')
load('PROMrun1.mat', 'fb')

atpfluxPROM = [vb(:,cell2mat(atpnocyclerxns(:,2)))]';
PROMatp = sum(bsxfun(@times,atpfluxPROM(1:end-1,:),cell2mat(atpnocyclerxns(1:end-1,3))));

soln = optimizeCbModel(iNJ661m);
FBAv = soln.x;
soln = optimizeCbModel(iNJ661m);
isequal(FBAv,soln.x);

FBAatp = [FBAv(cell2mat(atpnocyclerxns(:,2)))];
sum(cell2mat(atpnocyclerxns(1:end-1,3)).*FBAatp(1:end-1,:))

setdiff(atprxns(:,1),atpnocyclerxns(:,1))

FBAatp2 = [FBAv(cell2mat(atprxns(:,2)))];
sum(cell2mat(atprxns(1:end-1,3)).*FBAatp2(1:end-1,:))

atpfluxPROM2 = [vb(:,cell2mat(atprxns(:,2)))]';
PROMatp2 = sum(bsxfun(@times,atpfluxPROM2(1:end-1,:),cell2mat(atprxns(1:end-1,3))));

nadphrxns = findRxnsFromMets(iNJ661m,'nadph[c]');
nadphrxns(:,2) = cellfun(@(x) strmatch(x,iNJ661m.rxns,'exact'),nadphrxns(:,1),'UniformOutput',false);
nadphrxns(:,3) = cellfun(@(x) full(iNJ661m.S(181,x)),nadphrxns(:,2),'UniformOutput',false);

nadphfluxF = [v1(cell2mat(nadphrxns(:,2))) v2(cell2mat(nadphrxns(:,2)))];

sum(cell2mat(nadphrxns(1:end-1,3)).*nadphfluxF(1:end-1,1))
sum(cell2mat(nadphrxns(1:end-1,3)).*nadphfluxF(1:end-1,2))

FBAnadph2 = [FBAv(cell2mat(nadphrxns(:,2)))];
sum(cell2mat(nadphrxns(1:end-1,3)).*FBAnadph2(1:end-1,:))

nadphfluxPROM2 = [vb(:,cell2mat(nadphrxns(:,2)))]';
PROMnadph2 = sum(bsxfun(@times,nadphfluxPROM2(1:end-1,:),cell2mat(nadphrxns(1:end-1,3))));