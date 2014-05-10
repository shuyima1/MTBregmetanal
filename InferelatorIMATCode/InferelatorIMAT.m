function [InferelatorKOexpressionData,InferelatorKOexp, InferelatorIMATmodel, InferelatorIMATrxns] = InferelatorIMAT(model,MetTFCoefMat,MetTRNexp,thresh,runIMATflag,epsilon)

if ~exist('epsilon','var') || isempty(epsilon)
    epsilon = 1;
end

% Isolate the Inferelator Coefficient Matrix
CoefMat = cell2mat(MetTFCoefMat(2:end,2:end));

% Isolate the list of TF regulators
Regulators = MetTFCoefMat(1,2:end)';

% Isolate the list of Metabolic Gene targets
Targets = MetTFCoefMat(2:end,1);

% Initiate the output variables
InferelatorIMATmodel = cell(size(Regulators,1),1);
InferelatorIMATrxns = cell(size(Regulators,1),1);
InferelatorKOexp = zeros(size(Targets,1),size(Regulators,1));

[ia1, ia1, ib1] = intersect(Targets,MetTRNexp(:,1),'stable');
[ia2, ia2, ib2] = intersect(Regulators,MetTRNexp(:,1),'stable');

% Calculate the pre-knockout baseline expression (currently set to median across samples)
BaseExpressionProfile = median(cell2mat(MetTRNexp(:,2:end)),2);

% Extract the expression for the Metabolic target genes and the TFs (some
% values will not be available). These will be in the same order as the
% genes in 'Targets','Regulators'
BaseMetProfile = nan(size(Targets,1),1);
BaseMetProfile(ia1) = BaseExpressionProfile(ib1,:);

BaseTFProfile = nan(size(Regulators,1),1);
BaseTFProfile(ia2) = BaseExpressionProfile(ib2,:);

for nTF = 1:size(Regulators,1)
   
    %affectedMetGens lists the row numbers of the metabolic genes that are
    %affected by the TF knockout (ie. target genes that have a nonzero
    %coefficient of the TF in question
   affectedMetGens = find(CoefMat(:,nTF) ~= 0);
   
    % by default, set the predicted KO expression profile to be the same as
    % the BaseMetProfile, then modify by TF
   predictedKOexp = BaseMetProfile;
   
    % iterate through each of the affected target genes
   for tar = 1:numel(affectedMetGens)
       
       % find all of the TFs that determine the affected target gene's
       % expression
       nzTF = find(CoefMat(affectedMetGens(tar),:) ~= 0);
       
       % initiate the target gene's expression to 0
       predictedKOexp(affectedMetGens(tar)) = 0;
       
       % iterate through each of the TFs that determine the target gene's
       % expression
       for param = 1:numel(nzTF)
           
           % extract the expression of each of the TFs that determine the
           % target gene's expression. For the TF to be knocked out, set the
           % expression value to zero.
           if nzTF(param) == nTF
               tfexp = 0;
           else
               tfexp = BaseTFProfile(nzTF(param));
           end
           predictedKOexp(affectedMetGens(tar)) = predictedKOexp(affectedMetGens(tar)) + CoefMat(affectedMetGens(tar),nzTF(param))*tfexp;
       
       end
   end
   
   InferelatorKOexpressionData(nTF,1).Locus = Targets;
   
   %Binarize the expression data. Set Metabolic genes that could not be
   %predicted (because of nan TF values) to 1.
   InferelatorKOexpressionData(nTF,1).Data = [predictedKOexp > quantile(predictedKOexp,thresh)];
   InferelatorKOexpressionData(nTF,1).Data(isnan(predictedKOexp)) = 1;
   
   InferelatorKOexp(:,nTF) = predictedKOexp;
   
   if isequal(upper(runIMATflag),'IMAT')
       % Run iMAT on the resulting model
       [IMATmodel, IMATrxns] = createTissueSpecificModelimatorphanSM(model,InferelatorKOexpressionData(nTF,1),1,1,[],'iMAT',epsilon);
       InferelatorIMATmodel{nTF} = IMATmodel;
       InferelatorIMATrxns{nTF} = IMATrxns;
   elseif isequal(upper(runIMATflag),'GIMME')
       [IMATmodel, IMATrxns] = createTissueSpecificModel(model,InferelatorKOexpressionData(nTF,1),1,1,[],'GIMME');
       InferelatorIMATmodel{nTF} = IMATmodel;
       InferelatorIMATrxns{nTF} = IMATrxns;
   end
end