printRxnFormula(Beste7H9,'R864')
Beste7H9noaa = changeRxnBounds(Beste7H9,'R864',0,'u');
Beste7H9acnoaa = changeRxnBounds(Beste7H9ac,'R864',0,'u');
Beste7H9glcnoaa = changeRxnBounds(Beste7H9glc,'R864',0,'u');
Beste7H9glycnoaa = changeRxnBounds(Beste7H9glyc,'R864',0,'u');
Beste7H9ncnoaa = changeRxnBounds(Beste7H9nc,'R864',0,'u');

sB7_na = optimizeCbModel(Beste7H9noaa);
sB7ac_na = optimizeCbModel(Beste7H9acnoaa);
sB7glc_na = optimizeCbModel(Beste7H9glcnoaa);
sB7glyc_na = optimizeCbModel(Beste7H9glycnoaa);
sB7nc_na = optimizeCbModel(Beste7H9ncnoaa);

[grRatioB7noaa grKOB7noaa tmp tmp tmp sgfluxB7noaa] = singleGeneDeletion(Beste7H9noaa);
[grRatioB7acnoaa grKOB7acnoaa tmp tmp tmp sgfluxB7acnoaa] = singleGeneDeletion(Beste7H9acnoaa);
[grRatioB7glcnoaa grKOB7glcnoaa tmp tmp tmp sgfluxB7glcnoaa] = singleGeneDeletion(Beste7H9glcnoaa);
[grRatioB7glycnoaa grKOB7glycnoaa tmp tmp tmp sgfluxB7glycnoaa] = singleGeneDeletion(Beste7H9glycnoaa);
[grRatioB7ncnoaa grKOB7ncnoaa tmp tmp tmp sgfluxB7ncnoaa] = singleGeneDeletion(Beste7H9ncnoaa);

[fB7_na fkoB7_na vB7_na VkoB7_na sB7_na lostxnsB7_na probtfgeneB7_na] = promv2(Beste7H9noaa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7ac_na fkoB7ac_na vB7ac_na VkoB7ac_na sB7ac_na lostxnsB7ac_na probtfgeneB7ac_na] = promv2(Beste7H9acnoaa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7glc_na fkoB7glc_na vB7glc_na VkoB7glc_na sB7glc_na lostxnsB7glc_na probtfgeneB7glc_na] = promv2(Beste7H9glcnoaa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7glyc_na fkoB7glyc_na vB7glyc_na VkoB7glyc_na sB7glyc_na lostxnsB7glyc_na probtfgeneB7glyc_na] = promv2(Beste7H9glycnoaa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7nc_na fkoB7nc_na vB7nc_na VkoB7nc_na sB7nc_na lostxnsB7nc_na probtfgeneB7nc_na] = promv2(Beste7H9ncnoaa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
fB7_na = fB7_na';fB7ac_na = fB7ac_na';fB7glc_na = fB7glc_na';fB7glyc_na = fB7glyc_na';fB7nc_na = fB7nc_na';

load('BesteMediaModels.mat', 'BesteExchange')
aaXR = BesteExchange(20:23,1);
aaXR = [aaXR;BesteExchange([26;28],1)];
aaXR = [aaXR;BesteExchange([29:35],1)];
aaXR = [aaXR;BesteExchange([39 42 43 48 49 50 53 64 66 68:71 74 75 77 78 81 83 84 86 87 91:94 97 98],1)];
printRxnFormula(Beste7H9,aaXR)
aaXR(:,2) = ans;

Beste7H9aa = changeRxnBounds(Beste7H9,aaXR(:,1),1,'u');
Beste7H9acaa = changeRxnBounds(Beste7H9ac,aaXR(:,1),1,'u');
Beste7H9glcaa = changeRxnBounds(Beste7H9glc,aaXR(:,1),1,'u');
Beste7H9glycaa = changeRxnBounds(Beste7H9glyc,aaXR(:,1),1,'u');
Beste7H9ncaa = changeRxnBounds(Beste7H9nc,aaXR(:,1),1,'u');

sB7ac_aa = optimizeCbModel(Beste7H9acaa);
sB7glc_aa = optimizeCbModel(Beste7H9glcaa);
sB7glyc_aa = optimizeCbModel(Beste7H9glycaa);
sB7nc_aa = optimizeCbModel(Beste7H9ncaa);
sB7_aa = optimizeCbModel(Beste7H9aa);

[grRatioB7aa grKOB7aa tmp tmp tmp sgfluxB7aa] = singleGeneDeletion(Beste7H9aa);
[grRatioB7acaa grKOB7acaa tmp tmp tmp sgfluxB7acaa] = singleGeneDeletion(Beste7H9acaa);
[grRatioB7glcaa grKOB7glcaa tmp tmp tmp sgfluxB7glcaa] = singleGeneDeletion(Beste7H9glcaa);
[grRatioB7glycaa grKOB7glycaa tmp tmp tmp sgfluxB7glycaa] = singleGeneDeletion(Beste7H9glycaa);
[grRatioB7ncaa grKOB7ncaa tmp tmp tmp sgfluxB7ncaa] = singleGeneDeletion(Beste7H9ncaa);


[fB7_aa fkoB7_aa vB7_aa VkoB7_aa sB7_aa lostxnsB7_aa probtfgeneB7_aa] = promv2(Beste7H9aa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7ac_aa fkoB7ac_aa vB7ac_aa VkoB7ac_aa sB7ac_aa lostxnsB7ac_aa probtfgeneB7ac_aa] = promv2(Beste7H9acaa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7glc_aa fkoB7glc_aa vB7glc_aa VkoB7glc_aa sB7glc_aa lostxnsB7glc_aa probtfgeneB7glc_aa] = promv2(Beste7H9glcaa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7glyc_aa fkoB7glyc_aa vB7glyc_aa VkoB7glyc_aa sB7glyc_aa lostxnsB7glyc_aa probtfgeneB7glyc_aa] = promv2(Beste7H9glycaa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7nc_aa fkoB7nc_aa vB7nc_aa VkoB7nc_aa sB7nc_aa lostxnsB7nc_aa probtfgeneB7nc_aa] = promv2(Beste7H9ncaa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
fB7_aa = fB7_aa';fB7ac_aa = fB7ac_aa';fB7glc_aa = fB7glc_aa';fB7glyc_aa = fB7glyc_aa';fB7nc_aa = fB7nc_aa';
