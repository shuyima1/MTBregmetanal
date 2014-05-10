% BIGoperonints = {};
% for i = 1:size(MTBoperongens,1)
%     if size(MTBoperongens{i},2) > 1
%         for j = 1:size(MTBoperongens{i},2)
%             tmp = find(strcmp(MTBoperongens{i}{j},BIGtf_05_2013_Rv2621(:,2)));
%             for k = 1:size(tmp,1)
%                 BIGoperonints = [BIGoperonints; repmat(BIGtf_05_2013_Rv2621(tmp(k),1),size(MTBoperongens{i},2)-1,1) setdiff(MTBoperongens{i},MTBoperongens{i}{j})'];
%             end
%         end
%     end
% end
% BIGoperonints(:,2) = regexprep(BIGoperonints(:,2),'\s','');
%% Carbon Source
[f fko v vko stat lostxns probtfgene] = promv2activatoronly(Beste7H9modJ,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
B7JBIG50Activatorspromstats.f = f;
B7JBIG50Activatorspromstats.fko = fko;
B7JBIG50Activatorspromstats.v = v;
B7JBIG50Activatorspromstats.vko = vko;
B7JBIG50Activatorspromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2activatoronly(Beste7H9Jaa2,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
B7Jaa2BIG50Activatorspromstats.f = f;
B7Jaa2BIG50Activatorspromstats.fko = fko;
B7Jaa2BIG50Activatorspromstats.v = v;
B7Jaa2BIG50Activatorspromstats.vko = vko;
B7Jaa2BIG50Activatorspromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2activatoronly(BesteGriffinJ,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
BJGBIG50Activatorspromstats.f = f;
BJGBIG50Activatorspromstats.fko = fko;
BJGBIG50Activatorspromstats.v = v;
BJGBIG50Activatorspromstats.vko = vko;
BJGBIG50Activatorspromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2activatoronly(BesteSautonsJ,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
BJSBIG50Activatorspromstats.f = f;
BJSBIG50Activatorspromstats.fko = fko;
BJSBIG50Activatorspromstats.v = v;
BJSBIG50Activatorspromstats.vko = vko;
BJSBIG50Activatorspromstats.probtfgene = probtfgene;

% 7H9 acetate, glucose, glycerol
[f fko v vko stat lostxns probtfgene] = promv2activatoronly(Beste7H9Jac,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
B7JacBIG50Activatorspromstats.f = f;
B7JacBIG50Activatorspromstats.fko = fko;
B7JacBIG50Activatorspromstats.v = v;
B7JacBIG50Activatorspromstats.vko = vko;
B7JacBIG50Activatorspromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2activatoronly(Beste7H9Jglc,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
B7JglcBIG50Activatorspromstats.f = f;
B7JglcBIG50Activatorspromstats.fko = fko;
B7JglcBIG50Activatorspromstats.v = v;
B7JglcBIG50Activatorspromstats.vko = vko;
B7JglcBIG50Activatorspromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2activatoronly(Beste7H9Jglyc,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
B7JglycBIG50Activatorspromstats.f = f;
B7JglycBIG50Activatorspromstats.fko = fko;
B7JglycBIG50Activatorspromstats.v = v;
B7JglycBIG50Activatorspromstats.vko = vko;
B7JglycBIG50Activatorspromstats.probtfgene = probtfgene;

%sauton's acetate and glucose
[f fko v vko stat lostxns probtfgene] = promv2activatoronly(Beste7H9JSac,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
BJSacBIG50Activatorspromstats.f = f;
BJSacBIG50Activatorspromstats.fko = fko;
BJSacBIG50Activatorspromstats.v = v;
BJSacBIG50Activatorspromstats.vko = vko;
BJSacBIG50Activatorspromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2activatoronly(Beste7H9JSglc,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
BJSglcBIG50Activatorspromstats.f = f;
BJSglcBIG50Activatorspromstats.fko = fko;
BJSglcBIG50Activatorspromstats.v = v;
BJSglcBIG50Activatorspromstats.vko = vko;
BJSglcBIG50Activatorspromstats.probtfgene = probtfgene;
save PROMmodJBIGBIG50Activatorsanalresults *promstats
clear *promstats

%% aaModels
aaModelsSautonsBIG50Activators = cell(20,1);
for i = 1:20;
    aaModelsSautonsBIG50Activators{i,1} = changeRxnBounds(BesteSautonsJ,aaXR{i,1},1,'u');
%     aaModels{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
end

for i = 1:20
    %aaModels{i,2} = optimizeCbModel(aaModels{i,1});
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(aaModelsSautonsBIG50Activators{i,1});
    aaModelsSautonsBIG50Activators{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2activatoronly(aaModelsSautonsBIG50Activators{i,1},metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
    aaModelsSautonsBIG50Activators{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

aaJSBIG50Activatorsgrowthrateabs = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),20);
for i = 1:size(aaJSBIG50Activatorsgrowthrateabs,2)
    aaJSBIG50Activatorsgrowthrateabs(:,i) = aaModelsSautonsBIG50Activators{i,3}{1};
end

aaJSBIG50Activatorsgrowthratio = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),20);
for i = 1:size(aaJSBIG50Activatorsgrowthratio,2)
    aaJSBIG50Activatorsgrowthratio(:,i) = aaModelsSautonsBIG50Activators{i,3}{1}/aaModelsSautonsBIG50Activators{i,2}{2};
end
save aaModelsmodJSBIGBIG50Activators aaModelsSautonsBIG50Activators aaJSBIG50Activatorsgrowthrateabs aaJSBIG50Activatorsgrowthratio
clear aaModelsSautonsBIG50Activators

%%
aaModels7H9BIG50Activators = cell(20,1);
for i = 1:20;
    aaModels7H9BIG50Activators{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
%     aaModels{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
end

for i = 1:20
    %aaModels{i,2} = optimizeCbModel(aaModels{i,1});
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(aaModels7H9BIG50Activators{i,1});
    aaModels7H9BIG50Activators{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2activatoronly(aaModels7H9BIG50Activators{i,1},metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
    aaModels7H9BIG50Activators{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

aaJ7BIG50Activatorsgrowthrateabs = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),20);
for i = 1:size(aaJ7BIG50Activatorsgrowthrateabs,2)
    aaJ7BIG50Activatorsgrowthrateabs(:,i) = aaModels7H9BIG50Activators{i,3}{1};
end

aaJ7BIG50Activatorsgrowthratio = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),20);
for i = 1:size(aaJ7BIG50Activatorsgrowthratio,2)
    aaJ7BIG50Activatorsgrowthratio(:,i) = aaModels7H9BIG50Activators{i,3}{1}/aaModels7H9BIG50Activators{i,2}{2};
end

save aaModelsmodJ7BIGBIG50Activators aaModels7H9BIG50Activators aaJ7BIG50Activatorsgrowthratio aaJ7BIG50Activatorsgrowthrateabs
clear aaModels7H9BIG50Activators

%% AuxModels
AuxModelsSautonsBIG50Activators = cell(size(AuxExchanges,1),1);
for i = 1:size(AuxExchanges,1);
    AuxModelsSautonsBIG50Activators{i,1} = changeRxnBounds(BesteSautonsJ,AuxExchanges{i,1},1,'u');
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(AuxModelsSautonsBIG50Activators{i,1});
    AuxModelsSautonsBIG50Activators{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2activatoronly(AuxModelsSautonsBIG50Activators{i,1},metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
    AuxModelsSautonsBIG50Activators{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

AuxJSBIG50ActivatorsgrRatio = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJSBIG50ActivatorsgrRatio(:,i) = AuxModelsSautonsBIG50Activators{i,3}{1}/AuxModelsSautonsBIG50Activators{i,2}{2};end
AuxJSBIG50ActivatorsgrKO = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJSBIG50ActivatorsgrKO(:,i) = AuxModelsSautonsBIG50Activators{i,3}{1};end
save AuxModelsSautonsBIGBIG50Activators AuxModelsSautonsBIG50Activators AuxExchanges AuxJSBIG50ActivatorsgrRatio AuxJSBIG50ActivatorsgrKO
clear AuxModelsSautonsBIG50Activators 

AuxModels7H9Jaa2BIG50Activators = cell(size(AuxExchanges,1),1);
for i = 1:size(AuxExchanges,1);
    AuxModels7H9Jaa2BIG50Activators{i,1} = changeRxnBounds(Beste7H9Jaa2,AuxExchanges{i,1},1,'u');
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(AuxModels7H9Jaa2BIG50Activators{i,1});
    AuxModels7H9Jaa2BIG50Activators{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2activatoronly(AuxModels7H9Jaa2BIG50Activators{i,1},metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],[],[],0);
    AuxModels7H9Jaa2BIG50Activators{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

AuxJ7BIG50ActivatorsgrRatio = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJ7BIG50ActivatorsgrRatio(:,i) = AuxModels7H9Jaa2BIG50Activators{i,3}{1}/AuxModels7H9Jaa2BIG50Activators{i,2}{2};end
AuxJ7BIG50ActivatorsgrKO = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJ7BIG50ActivatorsgrKO(:,i) = AuxModels7H9Jaa2BIG50Activators{i,3}{1};end
save AuxModelsJ7aa2BIGBIG50Activators AuxModels7H9Jaa2BIG50Activators AuxExchanges AuxJ7BIG50ActivatorsgrRatio AuxJ7BIG50ActivatorsgrKO
clear AuxModels7H9Jaa2BIG50Activators 
%%
% fdoubleBIGBIG50ActivatorsJaa2 = PROMdoubleKO(Beste7H9Jaa2,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0);
% fdoubleBIGBIG50ActivatorsJS = PROMdoubleKO(BesteSautonsJ,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0);
% save Beste7H9Jaa2chipnetDoubleKOresults fdouble*
%clear fdoublechipnetJaa2
clear f fko v vko stat lostxns grRatio grWT grWT delRs delRs sgdflux