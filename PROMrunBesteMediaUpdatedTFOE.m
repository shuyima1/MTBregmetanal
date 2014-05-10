load('TFOE052813all.mat')
RvGeneIDs = geneIDs(9063:13061);
RvTFoeData = tfoeData(9063:13061,:);
save RvTFOE05283data RvGeneIDs RvTFoeData sampleIDs

initCobraToolbox

BIGtf_05_2013_Rv2621 = BIGtf_05_2013;
BIGtf_05_2013_Rv2621(x) = {'Rv2621'};
save mtbPROMmodelInputs BIGtf_05_2013_Rv2621 -append

MetGens = intersect(Beste7H9.genes,BIGtf_05_2013(:,2));
MetGens(:,2) = cellfun(@(x) strmatch(x,BIGtf_05_2013(:,2),'exact'),MetGens(:,1),'UniformOutput',false);
MetTFs = {};for i = 1:size(MetGens,1);for j = 1:size(MetGens{i,2},1);MetTFs = [MetTFs;BIGtf_05_2013(MetGens{i,2}(j),1)];end;end;MetTFs = unique(MetTFs);
MetSamples = cellfun(@(x) strmatch(x,sampleIDs(:,1)),MetTFs,'UniformOutput',false);
MetSamples{27,1} = strmatch('Rv2621',sampleIDs(:,1));

MetSamplesLong = [];for i = 1:44;MetSamplesLong = [MetSamplesLong;MetSamples{i,1}];end
MetSampleIDs = sampleIDs([[8:15]';MetSamplesLong],:);

metRvFC2 = tfoeFC2rv(:,MetSamplesLong);

[fB72e,f_koB72e,vB72e,v_koB72e,stat1B72e,lostxnsB72e,probtfgeneB72e] =  promv2(Beste7H9,metRvFC2,RvGeneIDs,BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],1);
[fB7ac2e,f_koB7ac2e,vB7ac2e,v_koB7ac2e,stat1B7ac2e,lostxnsB7ac2e,probtfgeneB7ac2e] =  promv2(Beste7H9ac,metRvFC2,RvGeneIDs,BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],1);
[fB7glc2e,f_koB7glc2e,vB7glc2e,v_koB7glc2e,status1B7glc2e,lostxnsB7glc2e,probtfgeneB7glc2e] =  promv2(Beste7H9glc,metRvFC2,RvGeneIDs,BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],1);
[fB7glyc2e,f_koB7glyc2e,vB7glyc2e,v_koB7glyc2e,status1B7glyc2e,lostxnsB7glyc2e,probtfgeneB7glyc2e] =  promv2(Beste7H9glyc,metRvFC2,RvGeneIDs,BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],1);
[fB7nc2e,f_koB7nc2e,vB7nc2e,v_koB7nc2e,status1B7nc2e,lostxnsB7nc2e,probtfgeneB7nc2e] =  promv2(Beste7H9nc,metRvFC2,RvGeneIDs,BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],1);
fB72e = fB72e';fB7ac2e = fB7ac2e';fB7glc2e = fB7glc2e';fB7glyc2e = fB7glyc2e';fB7nc2e = fB7nc2e';

sB7 = optimizeCbModel(Beste7H9,[],'one');
sB7ac = optimizeCbModel(Beste7H9ac,[],'one');
sB7glc = optimizeCbModel(Beste7H9glc,[],'one');
sB7glyc = optimizeCbModel(Beste7H9glyc,[],'one');
sB7nc = optimizeCbModel(Beste7H9nc,[],'one');

[grRatioB7 grKOB7 tmp tmp tmp sgfluxB7] = singleGeneDeletion(Beste7H9);
[grRatioB7ac grKOB7ac tmp tmp tmp sgfluxB7ac] = singleGeneDeletion(Beste7H9ac);
[grRatioB7glc grKOB7glc tmp tmp tmp sgfluxB7glc] = singleGeneDeletion(Beste7H9glc);
[grRatioB7glyc grKOB7glyc tmp tmp tmp sgfluxB7glyc] = singleGeneDeletion(Beste7H9glyc);
[grRatioB7nc grKOB7nc tmp tmp tmp sgfluxB7nc] = singleGeneDeletion(Beste7H9nc);

save PROMmediaBesteUpdatedTFOE120513 *B72e *B7ac2e *B7glc2e *B7glyc2e *B7nc2e MetSamplesLong metRvFC2 RvGeneIDs MetSampleIDs