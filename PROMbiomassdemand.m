bioMetTFess = cell(size(bioMets,1),2);
for i = 1:length(rxnNames)
    
    mB = changeObjective(modelBioDemand,rxnNames{i})
    [f,f_ko,v,v_ko,status1,lostxns,probtfgene] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);
    bioMetTFess{i,1} = f;
    bioMetTFess{i,2} = v;
end

save bioMetGenEss bioMetTFess