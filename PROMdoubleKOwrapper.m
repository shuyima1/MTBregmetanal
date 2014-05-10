fko = zeros(size(nonlethal,1),50);

for i = 1:size(nonlethal,1)
    [modelDel] = deleteModelGenes(Beste7H9,Beste7H9.genes{nonlethal(i)});
    [f] =  promv2(modelDel,metRvFC2,RvGeneIDs,BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],vmin,vmax,[],0,[],0);
    fko(i,:) = f;
end
