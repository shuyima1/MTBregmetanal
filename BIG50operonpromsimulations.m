BIGoperonints = {};
for i = 1:size(MTBoperongens,1)
    if size(MTBoperongens{i},2) > 1
        for j = 1:size(MTBoperongens{i},2)
            tmp = find(strcmp(MTBoperongens{i}{j},BIGtf_05_2013_Rv2621(:,2)));
            for k = 1:size(tmp,1)
                BIGoperonints = [BIGoperonints; repmat(BIGtf_05_2013_Rv2621(tmp(k),1),size(MTBoperongens{i},2)-1,1) setdiff(MTBoperongens{i},MTBoperongens{i}{j})'];
            end
        end
    end
end
BIGoperonints(:,2) = regexprep(BIGoperonints(:,2),'\s','');
%% Carbon Source
[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9modJ,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
B7J50OPpromstats.f = f;
B7J50OPpromstats.fko = fko;
B7J50OPpromstats.v = v;
B7J50OPpromstats.vko = vko;
B7J50OPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9Jaa2,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
B7Jaa250OPpromstats.f = f;
B7Jaa250OPpromstats.fko = fko;
B7Jaa250OPpromstats.v = v;
B7Jaa250OPpromstats.vko = vko;
B7Jaa250OPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(BesteGriffinJ,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
BJG50OPpromstats.f = f;
BJG50OPpromstats.fko = fko;
BJG50OPpromstats.v = v;
BJG50OPpromstats.vko = vko;
BJG50OPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(BesteSautonsJ,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
BJS50OPpromstats.f = f;
BJS50OPpromstats.fko = fko;
BJS50OPpromstats.v = v;
BJS50OPpromstats.vko = vko;
BJS50OPpromstats.probtfgene = probtfgene;

% 7H9 acetate, glucose, glycerol
[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9Jac,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
B7Jac50OPpromstats.f = f;
B7Jac50OPpromstats.fko = fko;
B7Jac50OPpromstats.v = v;
B7Jac50OPpromstats.vko = vko;
B7Jac50OPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9Jglc,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
B7Jglc50OPpromstats.f = f;
B7Jglc50OPpromstats.fko = fko;
B7Jglc50OPpromstats.v = v;
B7Jglc50OPpromstats.vko = vko;
B7Jglc50OPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9Jglyc,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
B7Jglyc50OPpromstats.f = f;
B7Jglyc50OPpromstats.fko = fko;
B7Jglyc50OPpromstats.v = v;
B7Jglyc50OPpromstats.vko = vko;
B7Jglyc50OPpromstats.probtfgene = probtfgene;

%sauton's acetate and glucose
[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9JSac,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
BJSac50OPpromstats.f = f;
BJSac50OPpromstats.fko = fko;
BJSac50OPpromstats.v = v;
BJSac50OPpromstats.vko = vko;
BJSac50OPpromstats.probtfgene = probtfgene;

[f fko v fko stat lostxns probtfgene] = promv2(Beste7H9JSglc,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
BJSglc50OPpromstats.f = f;
BJSglc50OPpromstats.fko = fko;
BJSglc50OPpromstats.v = v;
BJSglc50OPpromstats.vko = vko;
BJSglc50OPpromstats.probtfgene = probtfgene;
save PROMmodJBIG50OPanalresults *promstats
clear *promstats

%% aaModels
aaModelsSautons50OP = cell(20,1);
for i = 1:20;
    aaModelsSautons50OP{i,1} = changeRxnBounds(BesteSautonsJ,aaXR{i,1},1,'u');
%     aaModels{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
end

for i = 1:20
    %aaModels{i,2} = optimizeCbModel(aaModels{i,1});
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(aaModelsSautons50OP{i,1});
    aaModelsSautons50OP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(aaModelsSautons50OP{i,1},metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
    aaModelsSautons50OP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

aaJS50OPgrowthrateabs = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),20);
for i = 1:size(aaJS50OPgrowthrateabs,2)
    aaJS50OPgrowthrateabs(:,i) = aaModelsSautons50OP{i,3}{1};
end

aaJS50tOPgrowthratio = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),20);
for i = 1:size(aaJS50tOPgrowthratio,2)
    aaJS50tOPgrowthratio(:,i) = aaModelsSautons50OP{i,3}{1}/aaModelsSautons50OP{i,2}{2};
end
save aaModelsmodJSBIG50OP aaModelsSautons50OP aaJS50OPgrowthrateabs aaJS50tOPgrowthratio
clear aaModelsSautons50OP

%%
aaModels7H950OP = cell(20,1);
for i = 1:20;
    aaModels7H950OP{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
%     aaModels{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
end

for i = 1:20
    %aaModels{i,2} = optimizeCbModel(aaModels{i,1});
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(aaModels7H950OP{i,1});
    aaModels7H950OP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(aaModels7H950OP{i,1},metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
    aaModels7H950OP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

aaJ750OPgrowthrateabs = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),20);
for i = 1:size(aaJ750OPgrowthrateabs,2)
    aaJ750OPgrowthrateabs(:,i) = aaModels7H950OP{i,3}{1};
end

aaJ750OPgrowthratio = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),20);
for i = 1:size(aaJ750OPgrowthratio,2)
    aaJ750OPgrowthratio(:,i) = aaModels7H950OP{i,3}{1}/aaModels7H950OP{i,2}{2};
end

save aaModelsmodJ7BIG50OP aaModels7H950OP aaJ750OPgrowthratio aaJ750OPgrowthrateabs
clear aaModels7H950OP

%% AuxModels
AuxModelsSautons50OP = cell(size(AuxExchanges,1),1);
for i = 1:size(AuxExchanges,1);
    AuxModelsSautons50OP{i,1} = changeRxnBounds(BesteSautonsJ,AuxExchanges{i,1},1,'u');
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(AuxModelsSautons50OP{i,1});
    AuxModelsSautons50OP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(AuxModelsSautons50OP{i,1},metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
    AuxModelsSautons50OP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

AuxJS50OPgrRatio = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJS50OPgrRatio(:,i) = AuxModelsSautons50OP{i,3}{1}/AuxModelsSautons50OP{i,2}{2};end
AuxJS50OPgrKO = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJS50OPgrKO(:,i) = AuxModelsSautons50OP{i,3}{1};end
save AuxModelsSautonsBIG50OP AuxModelsSautons50OP AuxExchanges AuxJS50OPgrRatio AuxJS50OPgrKO
clear AuxModelsSautons50OP 

AuxModels7H9Jaa250OP = cell(size(AuxExchanges,1),1);
for i = 1:size(AuxExchanges,1);
    AuxModels7H9Jaa250OP{i,1} = changeRxnBounds(Beste7H9Jaa2,AuxExchanges{i,1},1,'u');
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(AuxModels7H9Jaa250OP{i,1});
    AuxModels7H9Jaa250OP{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(AuxModels7H9Jaa250OP{i,1},metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],[],[],0,[],0);
    AuxModels7H9Jaa250OP{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

AuxJ750OPgrRatio = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJ750OPgrRatio(:,i) = AuxModels7H9Jaa250OP{i,3}{1}/AuxModels7H9Jaa250OP{i,2}{2};end
AuxJ750OPgrKO = zeros(size(unique(BIGtf_05_2013_Rv2621(:,1)),1),size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJ750OPgrKO(:,i) = AuxModels7H9Jaa250OP{i,3}{1};end
save AuxModelsJ7aa2BIG50OP AuxModels7H9Jaa250OP AuxExchanges AuxJ750OPgrRatio AuxJ750OPgrKO
clear AuxModels7H9Jaa250OP 
%%
fdoubleBIG50OPJaa2 = PROMdoubleKO(Beste7H9Jaa2,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],0,[],0);
fdoubleBIG50OPJS = PROMdoubleKO(BesteSautonsJ,metRvFC2,RvGeneIDs(:,1),[BIGtf_05_2013_Rv2621(:,1);BIGoperonints(:,1)],[BIGtf_05_2013_Rv2621(:,2);BIGoperonints(:,2)],[],[],[],[],0,[],0);
save Beste7H9Jaa2chipnetDoubleKOresults fdouble*
%clear fdoublechipnetJaa2
clear f fko v vko stat lostxns grRatio grWT grWT delRs delRs sgdflux