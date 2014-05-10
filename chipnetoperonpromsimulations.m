% operonints = {};
% for i = 1:size(MTBoperongens,1)
%     if size(MTBoperongens{i},2) > 1
%         for j = 1:size(MTBoperongens{i},2)
%             tmp = find(strcmp(MTBoperongens{i}{j},chipnet(:,2)));
%             for k = 1:size(tmp,1)
%                operonints = [operonints; repmat(chipnet(tmp(k),1),size(MTBoperongens{i},2)-1,1) setdiff(MTBoperongens{i},MTBoperongens{i}{j})'];
%             end
%         end
%     end
% end
% operonints(:,2) = regexprep(operonints(:,2),'\s','');
%% Carbon Source
[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9modJ,RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
B7JcnetOPpromstats.f = f;
B7JcnetOPpromstats.fko = fko;
B7JcnetOPpromstats.v = v;
B7JcnetOPpromstats.vko = vko;
B7JcnetOPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9Jaa2,RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
B7Jaa2cnetOPpromstats.f = f;
B7Jaa2cnetOPpromstats.fko = fko;
B7Jaa2cnetOPpromstats.v = v;
B7Jaa2cnetOPpromstats.vko = vko;
B7Jaa2cnetOPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(BesteGriffinJ,RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
BJGcnetOPpromstats.f = f;
BJGcnetOPpromstats.fko = fko;
BJGcnetOPpromstats.v = v;
BJGcnetOPpromstats.vko = vko;
BJGcnetOPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(BesteSautonsJ,RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
BJScnetOPpromstats.f = f;
BJScnetOPpromstats.fko = fko;
BJScnetOPpromstats.v = v;
BJScnetOPpromstats.vko = vko;
BJScnetOPpromstats.probtfgene = probtfgene;

% 7H9 acetate, glucose, glycerol
[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9Jac,RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
B7JaccnetOPpromstats.f = f;
B7JaccnetOPpromstats.fko = fko;
B7JaccnetOPpromstats.v = v;
B7JaccnetOPpromstats.vko = vko;
B7JaccnetOPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9Jglc,RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
B7JglccnetOPpromstats.f = f;
B7JglccnetOPpromstats.fko = fko;
B7JglccnetOPpromstats.v = v;
B7JglccnetOPpromstats.vko = vko;
B7JglccnetOPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9Jglyc,RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
B7JglyccnetOPpromstats.f = f;
B7JglyccnetOPpromstats.fko = fko;
B7JglyccnetOPpromstats.v = v;
B7JglyccnetOPpromstats.vko = vko;
B7JglyccnetOPpromstats.probtfgene = probtfgene;

%sauton's acetate and glucose
[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9JSac,RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
BJSaccnetOPpromstats.f = f;
BJSaccnetOPpromstats.fko = fko;
BJSaccnetOPpromstats.v = v;
BJSaccnetOPpromstats.vko = vko;
BJSaccnetOPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9JSglc,RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
BJSglccnetOPpromstats.f = f;
BJSglccnetOPpromstats.fko = fko;
BJSglccnetOPpromstats.v = v;
BJSglccnetOPpromstats.vko = vko;
BJSglccnetOPpromstats.probtfgene = probtfgene;
save PROMmodJchipnetOPanalresults *promstats
clear *promstats

%% aaModels
aaModelsSautonschipnetOP = cell(20,1);
for i = 1:20;
    aaModelsSautonschipnetOP{i,1} = changeRxnBounds(BesteSautonsJ,aaXR{i,1},1,'u');
%     aaModels{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
end

for i = 1:20
    %aaModels{i,2} = optimizeCbModel(aaModels{i,1});
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(aaModelsSautonschipnetOP{i,1});
    aaModelsSautonschipnetOP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(aaModelsSautonschipnetOP{i,1},RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
    aaModelsSautonschipnetOP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

aaJScnetOPgrowthrateabs = zeros(size(unique(chipnet(:,1)),1),20);
for i = 1:size(aaJScnetOPgrowthrateabs,2)
    aaJScnetOPgrowthrateabs(:,i) = aaModelsSautonschipnetOP{i,3}{1};
end

aaJScnetOPgrowthratio = zeros(size(unique(chipnet(:,1)),1),20);
for i = 1:size(aaJScnetOPgrowthratio,2)
    aaJScnetOPgrowthratio(:,i) = aaModelsSautonschipnetOP{i,3}{1}/aaModelsSautonschipnetOP{i,2}{2};
end
save aaModelsmodJSchipnetOP aaModelsSautonschipnetOP aaJScnetOPgrowthrateabs aaJScnetOPgrowthratio
clear aaModelsSautonschipnetOP

%%
aaModels7H9chipnetOP = cell(20,1);
for i = 1:20;
    aaModels7H9chipnetOP{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
%     aaModels{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
end

for i = 1:20
    %aaModels{i,2} = optimizeCbModel(aaModels{i,1});
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(aaModels7H9chipnetOP{i,1});
    aaModels7H9chipnetOP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(aaModels7H9chipnetOP{i,1},RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
    aaModels7H9chipnetOP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

aaJ7cnetOPgrowthrateabs = zeros(size(unique(chipnet(:,1)),1),20);
for i = 1:size(aaJ7cnetOPgrowthrateabs,2)
    aaJ7cnetOPgrowthrateabs(:,i) = aaModels7H9chipnetOP{i,3}{1};
end

aaJ7cnetOPgrowthratio = zeros(size(unique(chipnet(:,1)),1),20);
for i = 1:size(aaJ7cnetOPgrowthratio,2)
    aaJ7cnetOPgrowthratio(:,i) = aaModels7H9chipnetOP{i,3}{1}/aaModels7H9chipnetOP{i,2}{2};
end

save aaModelsmodJ7chipnetOP aaModels7H9chipnetOP aaJ7cnetOPgrowthratio aaJ7cnetOPgrowthrateabs
clear aaModels7H9chipnetOP

%% AuxModels
AuxModelsSautonschipnetOP = cell(size(AuxExchanges,1),1);
for i = 1:size(AuxExchanges,1);
    AuxModelsSautonschipnetOP{i,1} = changeRxnBounds(BesteSautonsJ,AuxExchanges{i,1},1,'u');
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(AuxModelsSautonschipnetOP{i,1});
    AuxModelsSautonschipnetOP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(AuxModelsSautonschipnetOP{i,1},RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
    AuxModelsSautonschipnetOP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

AuxJScnetOPgrRatio = zeros(size(unique(chipnet(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJScnetOPgrRatio(:,i) = AuxModelsSautonschipnetOP{i,3}{1}/AuxModelsSautonschipnetOP{i,2}{2};end
AuxJScnetOPgrKO = zeros(size(unique(chipnet(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJScnetOPgrKO(:,i) = AuxModelsSautonschipnetOP{i,3}{1};end
save AuxModelsSautonschipnetOP AuxModelsSautonschipnetOP AuxExchanges AuxJScnetOPgrRatio AuxJScnetOPgrKO
clear AuxModelsSautonschipnetOP

AuxModels7H9Jaa2chipnetOP = cell(size(AuxExchanges,1),1);
for i = 1:size(AuxExchanges,1);
    AuxModels7H9Jaa2chipnetOP{i,1} = changeRxnBounds(Beste7H9Jaa2,AuxExchanges{i,1},1,'u');
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(AuxModels7H9Jaa2chipnetOP{i,1});
    AuxModels7H9Jaa2chipnetOP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(AuxModels7H9Jaa2chipnetOP{i,1},RvFCexp,Rvgens(:,1),[chipnet(:,1);operonints(:,1)],[chipnet(:,2);operonints(:,2)],[],[],[],[],[],[],0,[],0);
    AuxModels7H9Jaa2chipnetOP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

AuxJ7cnetOPgrRatio = zeros(size(unique(chipnet(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJ7cnetOPgrRatio(:,i) = AuxModels7H9Jaa2chipnetOP{i,3}{1}/AuxModels7H9Jaa2chipnetOP{i,2}{2};end
AuxJ7cnetOPgrKO = zeros(size(unique(chipnet(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJ7cnetOPgrKO(:,i) = AuxModels7H9Jaa2chipnetOP{i,3}{1};end
save AuxModelsJaa2chipnetOP AuxModels7H9Jaa2chipnetOP AuxExchanges AuxJ7cnetOPgrRatio AuxJ7cnetOPgrKO
clear AuxModels7H9Jaa2chipnetOP
%%
%fdoublechipnetJaa2 = PROMdoubleKO(Beste7H9Jaa2,RvFCexp,Rvgens(:,1),chipnet(:,1),chipnet(:,2),[],[],[],[],0,[],0);
%save Beste7H9Jaa2chipnetDoubleKOresults fdoublechipnetJaa2
%clear fdoublechipnetJaa2
clear f fko v vko stat lostxns grRatio grWT grWT delRs delRs sgdflux