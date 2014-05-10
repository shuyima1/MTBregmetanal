function [fko] = PROMdoubleKO(model,expression,expressionid,regulator,targets,litevidence,prob_prior,subsets,KAPPA,datathresh,probtfgene,sizeflag)

numtf = size(unique(regulator),1);
fko = nan(size(model.genes,1),numtf);

[grRatio grKO] = singleGeneDeletion(model);
nonlethal = find(grRatio > 0.01);

fko(grRatio <= 0.01,:) = repmat(grKO(grRatio <= 0.01),1,numtf);
clear grKO

for i = 1:size(nonlethal,1)
    [modelDel] = deleteModelGenes(model,model.genes{nonlethal(i)});
    [f] =  promv2(modelDel,expression,expressionid,regulator,targets,litevidence,prob_prior,subsets,[],[],KAPPA,datathresh,probtfgene,sizeflag);
    fko(nonlethal(i),:) = f;
end
