initCobraToolbox
% load('BesteGriffin.mat')
%load('Beste7H9aa2.mat')
% load('BesteMediaModels.mat', 'BesteRxns')
%load('BesteMediaModels.mat', 'Beste7H9')
%load('GriffinData.mat')
%load('mtbPROMinputs2.mat', 'metRvFC2')
%load('mtbPROMinputs2.mat', 'BIGtf_05_2013_Rv2621')
%load('mtbPROMinputs2.mat', 'RvGeneIDs')
%load('tfoeindivarraysRvdata.mat')

%RvFC = bsxfun(@minus,RvExp(:,16:end),RvExp(:,1));
%samplesTFOE= samples(16:end);
%clear samples
%uTFOE = unique(samplesTFOE);

%[c ia ib] = intersect(RvGeneIDs,Beste7H9.genes);

%TFOEx.Locus = RvGeneIDs(ia);
%RvFCmet = RvFC(ia,:);
%clear RvFC

%load MTBtfoeIMATinputs
samplesRv = regexp(samplesTFOE,'(\-|\_)','split');
samplesRv = cellfun(@(x) x{1},samplesRv,'UniformOutput',false);
uTFOE = unique(samplesRv);
TFOEx.Locus = metRvGeneIDs;
TFOErxnResults = cell(size(uTFOE));

for tf = 1:size(uTFOE,1)
    datatmp = RvFCmet(:,strmatch(uTFOE{tf},samplesTFOE));
    databin0 = datatmp < 0;
    databin0 = sum(databin0,2)./size(databin0,2) > 0.75;
    TFOEx.Data = ~databin0;
    try
        [model rxns] = createTissueSpecificModelimatorphanSM(Beste7H9aa2,TFOEx,1,1,[],'iMAT',0.001);
        TFOErxnResults{tf,1} = rxns;
    catch err
        TFOErxnresults{tf,1} = nan;
    end
end

TFOErxnResults(:,2) = cellfun(@(x) x.solution.cont(724),TFOErxnResults(:,1),'UniformOutput',false);
TFOErxnResults(:,3) = uTFOE(:,1);

save MTBtfoeIMATresults TFOErxnResults