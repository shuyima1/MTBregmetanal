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

%load AltPeakNetwork AltPeakOperonNet
%load Beste7H9modJ
load PROMchipnetInputs RvFCexp Rvgens AuxExchanges


%initCobraToolbox
%changeCobraSolver('glpk')

%% Carbon Source
[f fko v vko stat lostxns probtfgene] = promv2(Beste7H9modJ,RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
B7JAltPeakOPpromstats.f = f;
B7JAltPeakOPpromstats.fko = fko;
B7JAltPeakOPpromstats.v = v;
B7JAltPeakOPpromstats.vko = vko;
B7JAltPeakOPpromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2(Beste7H9Jaa2,RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
B7Jaa2AltPeakOPpromstats.f = f;
B7Jaa2AltPeakOPpromstats.fko = fko;
B7Jaa2AltPeakOPpromstats.v = v;
B7Jaa2AltPeakOPpromstats.vko = vko;
B7Jaa2AltPeakOPpromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2(BesteGriffinJ,RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
BJGAltPeakOPpromstats.f = f;
BJGAltPeakOPpromstats.fko = fko;
BJGAltPeakOPpromstats.v = v;
BJGAltPeakOPpromstats.vko = vko;
BJGAltPeakOPpromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2(BesteSautonsJ,RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
BJSAltPeakOPpromstats.f = f;
BJSAltPeakOPpromstats.fko = fko;
BJSAltPeakOPpromstats.v = v;
BJSAltPeakOPpromstats.vko = vko;
BJSAltPeakOPpromstats.probtfgene = probtfgene;

% 7H9 acetate, glucose, glycerol
[f fko v vko stat lostxns probtfgene] = promv2(Beste7H9Jac,RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
B7JacAltPeakOPpromstats.f = f;
B7JacAltPeakOPpromstats.fko = fko;
B7JacAltPeakOPpromstats.v = v;
B7JacAltPeakOPpromstats.vko = vko;
B7JacAltPeakOPpromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2(Beste7H9Jglc,RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
B7JglcAltPeakOPpromstats.f = f;
B7JglcAltPeakOPpromstats.fko = fko;
B7JglcAltPeakOPpromstats.v = v;
B7JglcAltPeakOPpromstats.vko = vko;
B7JglcAltPeakOPpromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2(Beste7H9Jglyc,RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
B7JglycAltPeakOPpromstats.f = f;
B7JglycAltPeakOPpromstats.fko = fko;
B7JglycAltPeakOPpromstats.v = v;
B7JglycAltPeakOPpromstats.vko = vko;
B7JglycAltPeakOPpromstats.probtfgene = probtfgene;

%sauton's acetate and glucose
[f fko v vko stat lostxns probtfgene] = promv2(Beste7H9JSac,RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
BJSacAltPeakOPpromstats.f = f;
BJSacAltPeakOPpromstats.fko = fko;
BJSacAltPeakOPpromstats.v = v;
BJSacAltPeakOPpromstats.vko = vko;
BJSacAltPeakOPpromstats.probtfgene = probtfgene;

[f fko v vko stat lostxns probtfgene] = promv2(Beste7H9JSglc,RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
BJSglcAltPeakOPpromstats.f = f;
BJSglcAltPeakOPpromstats.fko = fko;
BJSglcAltPeakOPpromstats.v = v;
BJSglcAltPeakOPpromstats.vko = vko;
BJSglcAltPeakOPpromstats.probtfgene = probtfgene;
save PROMmodJAltPeakOPanalresults *promstats
clear *promstats

%% aaModels
aaModelsSautonsAltPeakOP = cell(20,1);
for i = 1:20;
    aaModelsSautonsAltPeakOP{i,1} = changeRxnBounds(BesteSautonsJ,aaXR{i,1},1,'u');
%     aaModels{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
end

for i = 1:20
    %aaModels{i,2} = optimizeCbModel(aaModels{i,1});
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(aaModelsSautonsAltPeakOP{i,1});
    aaModelsSautonsAltPeakOP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(aaModelsSautonsAltPeakOP{i,1},RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
    aaModelsSautonsAltPeakOP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

aaJSAltPeakOPgrowthrateabs = zeros(size(unique(AltPeakOperonNet(:,1)),1),20);
for i = 1:size(aaJSAltPeakOPgrowthrateabs,2)
    aaJSAltPeakOPgrowthrateabs(:,i) = aaModelsSautonsAltPeakOP{i,3}{1};
end

aaJSAltPeakOPgrowthratio = zeros(size(unique(AltPeakOperonNet(:,1)),1),20);
for i = 1:size(aaJSAltPeakOPgrowthratio,2)
    aaJSAltPeakOPgrowthratio(:,i) = aaModelsSautonsAltPeakOP{i,3}{1}/aaModelsSautonsAltPeakOP{i,2}{2};
end
save aaModelsmodJSAltPeakOP aaModelsSautonsAltPeakOP aaJSAltPeakOPgrowthrateabs aaJSAltPeakOPgrowthratio
clear aaModelsSautonsAltPeakOP

%%
aaModels7H9AltPeakOP = cell(20,1);
for i = 1:20;
    aaModels7H9AltPeakOP{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
%     aaModels{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
end

for i = 1:20
    %aaModels{i,2} = optimizeCbModel(aaModels{i,1});
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(aaModels7H9AltPeakOP{i,1});
    aaModels7H9AltPeakOP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(aaModels7H9AltPeakOP{i,1},RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
    aaModels7H9AltPeakOP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

aaJ7AltPeakOPgrowthrateabs = zeros(size(unique(AltPeakOperonNet(:,1)),1),20);
for i = 1:size(aaJ7AltPeakOPgrowthrateabs,2)
    aaJ7AltPeakOPgrowthrateabs(:,i) = aaModels7H9AltPeakOP{i,3}{1};
end

aaJ7AltPeakOPgrowthratio = zeros(size(unique(AltPeakOperonNet(:,1)),1),20);
for i = 1:size(aaJ7AltPeakOPgrowthratio,2)
    aaJ7AltPeakOPgrowthratio(:,i) = aaModels7H9AltPeakOP{i,3}{1}/aaModels7H9AltPeakOP{i,2}{2};
end

save aaModelsmodJ7AltPeakOP aaModels7H9AltPeakOP aaJ7AltPeakOPgrowthratio aaJ7AltPeakOPgrowthrateabs
clear aaModels7H9AltPeakOP

%% AuxModels
AuxModelsSautonsAltPeakOP = cell(size(AuxExchanges,1),1);
for i = 1:size(AuxExchanges,1);
    AuxModelsSautonsAltPeakOP{i,1} = changeRxnBounds(BesteSautonsJ,AuxExchanges{i,1},1,'u');
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(AuxModelsSautonsAltPeakOP{i,1});
    AuxModelsSautonsAltPeakOP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(AuxModelsSautonsAltPeakOP{i,1},RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
    AuxModelsSautonsAltPeakOP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

AuxJSAltPeakOPgrRatio = zeros(size(unique(AltPeakOperonNet(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJSAltPeakOPgrRatio(:,i) = AuxModelsSautonsAltPeakOP{i,3}{1}/AuxModelsSautonsAltPeakOP{i,2}{2};end
AuxJSAltPeakOPgrKO = zeros(size(unique(AltPeakOperonNet(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJSAltPeakOPgrKO(:,i) = AuxModelsSautonsAltPeakOP{i,3}{1};end
save AuxModelsSautonsAltPeakOP AuxModelsSautonsAltPeakOP AuxExchanges AuxJSAltPeakOPgrRatio AuxJSAltPeakOPgrKO
clear AuxModelsSautonsAltPeakOP

AuxModels7H9Jaa2AltPeakOP = cell(size(AuxExchanges,1),1);
for i = 1:size(AuxExchanges,1);
    AuxModels7H9Jaa2AltPeakOP{i,1} = changeRxnBounds(Beste7H9Jaa2,AuxExchanges{i,1},1,'u');
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(AuxModels7H9Jaa2AltPeakOP{i,1});
    AuxModels7H9Jaa2AltPeakOP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(AuxModels7H9Jaa2AltPeakOP{i,1},RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);
    AuxModels7H9Jaa2AltPeakOP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

AuxJ7AltPeakOPgrRatio = zeros(size(unique(AltPeakOperonNet(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJ7AltPeakOPgrRatio(:,i) = AuxModels7H9Jaa2AltPeakOP{i,3}{1}/AuxModels7H9Jaa2AltPeakOP{i,2}{2};end
AuxJ7AltPeakOPgrKO = zeros(size(unique(AltPeakOperonNet(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJ7AltPeakOPgrKO(:,i) = AuxModels7H9Jaa2AltPeakOP{i,3}{1};end
save AuxModelsJaa2AltPeakOP AuxModels7H9Jaa2AltPeakOP AuxExchanges AuxJ7AltPeakOPgrRatio AuxJ7AltPeakOPgrKO
clear AuxModels7H9Jaa2AltPeakOP
%%

ExtractTFinfoAltPeakOperonNetJsauton

%fdoubleAltPeakOJaa2 = PROMdoubleKO(Beste7H9Jaa2,RvFCexp,Rvgens(:,1),chipnet(:,1),chipnet(:,2),[],[],[],[],0,[],0);
%save Beste7H9Jaa2chipnetDoubleKOresults fdoublechipnetJaa2
%clear fdoublechipnetJaa2
clear f fko v vko stat lostxns grRatio grWT grWT delRs delRs sgdflux