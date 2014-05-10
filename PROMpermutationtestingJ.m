numperm = 500;

%[vmin7 vmax7] = fluxVariability(Beste7H9Jaa2);

fnullexppermJ = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),numperm);
for i = 1:numperm
    expperm = metRvFC2(randperm(size(metRvFC2,1)),:);
    fnullexppermJ(:,i) = promv2activatoronly(Beste7H9Jaa2,expperm,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],vmin7,vmax7,[],[],[],0);
    %fnullexppermG(:,i) = promv2(BesteGriffinJ,expperm,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],vminG,vmaxG,[],[],[],0);
end
save PROMB7Jaa2permutationactivator fnullexppermJ

fnullnetpermJ = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),numperm);
for i = 1:numperm
    netperm = [BIGtf_05_2013_Rv2621(:,1) BIGtf_05_2013_Rv2621(randperm(size(BIGtf_05_2013_Rv2621,1)),2)];
    fnullnetpermJ(:,i) = promv2activatoronly(Beste7H9Jaa2,metRvFC2,RvGeneIDs(:,1),netperm(:,1),netperm(:,2),[],[],[],vmin7,vmax7,[],[],[],0);
    %fnullnetpermG(:,i) = promv2(BesteGriffinJ,metRvFC2,RvGeneIDs(:,1),netperm(:,1),netperm(:,2),[],[],[],vminG,vmaxG,[],[],[],0);
end
save PROMB7Jaa2permutationactivator fnullnetpermJ -append


    
