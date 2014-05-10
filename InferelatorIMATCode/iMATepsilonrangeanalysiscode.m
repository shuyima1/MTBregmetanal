load('yeastShlomiIMATmedianexp.mat', 'yeast6')
load('YeastPROMinputs.mat', 'YeastMetTRNexp')

initCobraToolbox

[vmin vmax] = fluxVariability(yeast6);
x = find(vmin == vmax);
xpos = find(vmin > 0 & vmax > 0);
xneg = find(vmin < 0 & vmax < 0);
negflux = [xneg vmin(xneg) vmax(xneg)];
posflux = [xpos vmin(xpos) vmax(xpos)];

posfluxsort = sort(posflux(:,2),'ascend');
negfluxsort = sort(negflux(:,3),'descend');

data2.Locus = YeastMetTRNexp(:,1);
data2.Data = cell2mat(YeastMetTRNexp(:,2))>7;

strmatch('growth',yeast6.rxnNames) % = 1876
epsilonrange = [0.001:0.001:0.02]';
rxnsrange = cell(numel(epsilonrange),1);

for i = 1:numel(epsilonrange)
    [rxns, rxns] = createTissueSpecificModelimatorphanSM(yeast6,data2,1,1,[],'iMAT',epsilonrange(i));
    rxnsrange{i} = rxns;
end

epsgrowthrange = zeros(size(epsilonrange));
for i = 1:numel(epsilonrange)
    epsgrowthrange(i)=rxnsrange{i}.solution.cont(1876);
end

save yeastIMATresults epsilonrange rxnsrange epsgrowthrange

solny6 = optimizeCbModel(yeast6,[],'one')
load('YeastMetCoef.mat', 'MetTFCoefMat')
[InferelatorKOexpressionDataQeps InfKOexpQeps InferelatorIMATmodelQeps InferelatorIMATrxnsQeps] = InferelatorIMAT(yeast6,MetTFCoefMat,YeastMetTRNexp,0.15,'IMAT',0.003);
fluxQeps = zeros(size(yeast6.rxns,1),size(InferelatorIMATrxnsQeps,1));
for i = 1:size(InferelatorIMATrxnsQeps,1)
    fluxQeps(:,i) = InferelatorIMATrxnsQeps{i}.solution.cont;
end

growthQeps = fluxQeps(1876,:)';
growthratioQeps = growthQeps/epsgrowthrange(3);growthratioQeps = [MetTFCoefMat(1,2:end)' num2cell(growthratioQeps)];
save yeastInferelatorIMATresults *Qeps epsilonrange


InfIMATrxnsEPS = cell(size(epsilonrange));
for i = 1:size(epsilonrange,1)
    [rxns rxns] = expressionmatrixIMAT(yeast6,InferelatorKOexpressionDataQeps(1).Locus,InfKOexpQeps,0.15,'IMAT',epsilonrange(i));
    InfIMATrxnsEPS{i} = rxns;
end

InfIMATgrEPS = nan(91,20);
for i = 1:20
    for j = 1:91
        if ~isempty(InfIMATrxnsEPS{i}{j})
            InfIMATgrEPS(j,i) = InfIMATrxnsEPS{i}{j}.solution.cont(1876);
        end
    end
end

save yeastInferelatorIMATepsrangeresults *EPS epsilonrange

figure;imagesc(InfIMATgrEPS)
colormap('hot')