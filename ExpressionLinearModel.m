% load('tfoeindivarraysRvdata.mat', 'RvExp')
% load('tfoeindivarraysRvdata.mat', 'samples')
% load('Beste7H9modJ.mat', 'Beste7H9Jaa2')
% 
% uTF = unique(BIGtf_05_2013_Rv2621(:,1));
% uTar = unique(BIGtf_05_2013_Rv2621(:,2));
% 
% uMetTar = intersect(BIGtf_05_2013_Rv2621(:,2),Beste7H9Jaa2.genes);
% uMetTar(:,2) = cellfun(@(x) BIGtf_05_2013_Rv2621(strcmp(x,BIGtf_05_2013_Rv2621(:,2)),1),uMetTar(:,1),'UniformOutput',false);
% 
% [uMetTarIX uMetTarIX] = intersect(Rvgens,uMetTar(:,1),'stable');
% uMetTarIX = num2cell(uMetTarIX);
% 
% for i = 1:size(uMetTarIX,1);uMetTarIX{i,2} = cellfun(@(x) find(strcmp(x,Rvgens)),uMetTar{i,2});end
% save MTBtarexpressionmodelInputs uMetTar*

% yMetTarExp = RvExp(cell2mat(uMetTarIX(:,1)),:);
expmodels = cell(size(uMetTarIX,1),2);
for i = 1:size(uMetTarIX,1)
    nparam = numel(uMetTarIX{i,2});
    mdl = LinearModel.fit(RvExp(uMetTarIX{i,2},:)',RvExp(uMetTarIX{i,1},:)');
    expmodels{i,1} = mdl;
    modelstats1 = anova(mdl);
    goodparams = find(modelstats1.pValue(1:nparam) < 0.05);
    if numel(goodparams) > 1
        nk = nchoosek(goodparams,2);
        Tmat = zeros(numel(goodparams)+1,nparam);
        for j = 1:numel(goodparams)
            Tmat(j+1,goodparams(j)) = 1; 
        end

        Tmatcum = Tmat;

        for j = 1:size(nk,1)
            Tmat1 = zeros(1,nparam);
            Tmat1(nk(j,:)) = 1;
            mdl2 = LinearModel.fit(RvExp(uMetTarIX{i,2},:)',RvExp(uMetTarIX{i,1},:)',[Tmat;Tmat1]);
            intstats = anova(mdl2);
            if intstats.pValue(numel(goodparams)+1) < 0.05
                Tmatcum = [Tmatcum;Tmat1];
            end
        end

        expmodels{i,2} = LinearModel.fit(RvExp(uMetTarIX{i,2},:)',RvExp(uMetTarIX{i,1},:)',Tmatcum);
        modelstats2 = anova(expmodels{i,2});
        if size(Tmatcum,1) > numel(goodparams)+1
            badparams = find(modelstats2.pValue(numel(goodparams)+1:end-1) >= 0.05);
            if ~isempty(badparams)
                Tmatcum(size(Tmat,1)+badparams,:) = [];
                expmodels{i,2} = LinearModel.fit(RvExp(uMetTarIX{i,2},:)',RvExp(uMetTarIX{i,1},:)',Tmatcum);
            end
        end
    end
end

