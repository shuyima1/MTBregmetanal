initCobraToolbox

epsilon = 0.003; % insert the epsilon value here
grix = 1876; %insert the row number of the growth reaction here

%% Run IMAT on each of the TFKOs
TFKOexpbin = TFKOTargetGenesLogRatio < 0;
TFKOx.Locus = Genes;
TFKOimatresults = cell(size(TF269,1),2);
TFKOimatgrowthrate = zeros(size(TF269,1),1);
for gen = 1:size(TF269,1)
   TFKOx.Data = TFKOexpbin(:,gen); 
   [model results] = createTissueSpecificModelimatorphanSM(yeast6,TFKOx,1,1,[],'iMAT',epsilon);
   TFKOimatresults{gen,1} = results;
   TFKOimatresults{gen,2} = model;
   TFKOimatgrowthrate(gen) = results.solution.cont(grix);
end

%% Run InferelatorIMAT on each of the TFKOs


