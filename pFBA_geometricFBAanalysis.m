initCobraToolbox

load('FangCycleMetAnalysis.mat', 'atprxns')
load('FangCycleMetAnalysis.mat', 'atpnocyclerxns')
load('FangCycleMetAnalysis.mat', 'negrxns')
load('FangCycleMetAnalysis.mat', 'posrxns')

% isolating the atp-specific reactions
FBAv = soln.x;
FBAatp = [FBAv(cell2mat(atpnocyclerxns(:,2)))];
sum(cell2mat(atpnocyclerxns(1:end-1,3)).*FBAatp(1:end-1,:))
setdiff(atprxns(:,1),atpnocyclerxns(:,1))

FBAatp2 = [FBAv(cell2mat(atprxns(:,2)))];
sum(cell2mat(atprxns(1:end-1,3)).*FBAatp2(1:end-1,:))

sum(cell2mat(atpnocyclerxns(1:end,3)).*FBAatp(1:end,:))
sum(cell2mat(atprxns(1:end,3)).*FBAatp2(1:end,:))

% generating the pFBA model for Fang
iNJ661m_pFBAanalysis

% checking for unconstrained reactions 
exfm = strmatch('EX_',modelIrrevFM2.rxns);
modelIrrefFM2nx = changeRxnBounds(modelIrrevFM2,modelIrrevFM2.rxns(exfm),0,'b');
modelIrrefFM2nx = changeRxnBounds(modelIrrefFM2nx,modelIrrevFM2.rxns(1087),0,'l');
[v1fmnx v2fmnx] = fluxVariability(modelIrrefFM2nx);

% generating geometric FBA solutions for Fang (had trouble with error:
% non-convergence
gFBAv = geometricFBA(iNJ661m),
gFBAredv = geometricFBA(mreduced);
optimizeCbModel(mreduced)

% generating geometric FBA solutions for jamshidi
load('MTBmodels.mat')
gFBAjamshidiv = geometricFBA(jamshidiBiGG);
mrxnsJam = cellfun(@(x) [find(jamshidiBiGG.S(x,:) ~= 0)' full(jamshidiBiGG.S(x,jamshidiBiGG.S(x,:) ~= 0))'],num2cell([1:size(jamshidiBiGG.mets,1)]'),'UniformOutput',false);
mrxnsJam(:,2) = cellfun(@(x) gFBAjamshidiv(x(:,1)),mrxnsJam(:,1),'UniformOutput',false);
mrxnsJam(:,3) = cellfun(@(x) sum(gFBAjamshidiv(x(:,1)).*x(:,2)),mrxnsJam(:,1),'UniformOutput',false);
fluxsumJamgFBA = cell2mat(mrxnsJam(:,3));

x = find(abs(fluxsumJamgFBA) > 10^-5);
x = strmatch('atp[c]',jamshidiBiGG.mets);
x
sum(gFBAjamshidiv(mrxnsJam{180,1}(1:end-1,1)).*mrxnsJam{180,1}(1:end-1,2))

% comparing the minimizing taxi-cab distance setting in optimizeCbModel
% with the geometric FBA results
solnone = optimizeCbModel(jamshidiBiGG,[],'one');
isequal(solnone.x,gFBAjamshidiv)
hist(solnone.x-gFBAjamshidiv)
solnoneFang = optimizeCbModel(iNJ661m,[],'one');