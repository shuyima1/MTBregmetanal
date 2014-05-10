%%% This code examines the effect of PROM KOs on nad/nadh and nadp/nadph
%%% demand reaction fluxes from the hypoxia models
load('iNJ661mModel.mat', 'iNJ661m')
biox = find(iNJ661m.c == 1);
bioS = iNJ661m.S(:,biox);
bioSm = find(bioS ~= 0);
bioMets = iNJ661m.mets(bioSm);
[modelBioDemand rxnNames] = addDemandReaction(iNJ661m,bioMets);
bioMetGenEss = cell(size(bioMets));
for i = 1:length(rxnNames);mB = changeObjective(modelBioDemand,rxnNames{i});[grRatio grRateKO grRateWT] = singleGeneDeletion(mB);bioMetGenEss{i} = grRateKO;end
save bioMetGenEss bioMetGenEss bioMets modelBioDemand rxnNames

mB = changeObjective(H1modelDemand,H1rxnNames{65});[grRatioH1nad grRateKOH1nad grRateWTH1nad] = singleGeneDeletion(mB);[fH1nad,f_koH1nad,vH1nad,v_ko,status1H1nad,lostxnsH1nad,probtfgeneH1nad] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);
mB = changeObjective(H1modelDemand,H1rxnNames{66});[grRatioH1nadp grRateKOH1nadp grRateWTH1nadp] = singleGeneDeletion(mB);[fH1nadp,f_koH1nadp,vH1nadp,v_koH1nadp,status1H1nadp,lostxnsH1nadp,probtfgeneH1nadp] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);

strmatch('nadh',H1modelDemand.mets)
strmatch('nadph',H1modelDemand.mets)
[modelBioDemandnadh rxnNamesNADH] = addDemandReaction(H1modelDemand,H1modelDemand.mets([513 515]));
mB = changeObjective(modelBioDemandnadh,rxnNamesNADH{1});[grRatioH1nadh grRateKOH1nadh grRateWTH1nadh] = singleGeneDeletion(mB);[fH1nadh,f_koH1nadh,vH1nadh,v_koH1nadh,status1H1nadh,lostxnsH1nadh,probtfgeneH1nadh] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);
mB = changeObjective(modelBioDemandnadh,rxnNamesNADH{2});[grRatioH1nadph grRateKOH1nadph grRateWTH1nadph] = singleGeneDeletion(mB);[fH1nadph,f_koH1nadph,vH1nadph,v_koH1nadph,status1H1nadph,lostxnsH1nadph,probtfgeneH1nadph] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);
load('bioMetGenEss.mat', 'modelBioDemand')
H1modelBioDemandnadh = modelBioDemandnadh;
H1rxnNamesNADH = rxnNamesNADH;
[modelBioDemandnadh rxnNamesNADH] = addDemandReaction(modelBioDemand,modelBioDemand.mets([513 515]));

strmatch('nadph',modelBioDemand.mets)
strmatch('nadh',modelBioDemand.mets)
[modelBioDemandnadh rxnNamesNADH] = addDemandReaction(modelBioDemand,modelBioDemand.mets([536 538]));
mB = changeObjective(modelBioDemandnadh,rxnNamesNADH{1});[grRationadh grRateKOnadh grRateWTnadh] = singleGeneDeletion(mB);[fnadh,f_konadh,vnadh,v_konadh,status1nadh,lostxnsnadh,probtfgenenadh] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);
mB = changeObjective(modelBioDemandnadh,rxnNamesNADH{2});[grRationadph grRateKOnadph grRateWTnadph] = singleGeneDeletion(mB);[fnadph,f_konadph,vnadph,v_konadph,status1nadph,lostxnsnadph,probtfgenenadph] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);

strmatch('nadh',H7modelDemand.mets)
strmatch('nadph',H7modelDemand.mets)
[H7modelBioDemandnadh H7rxnNamesNADH] = addDemandReaction(H7modelDemand,H7modelDemand.mets([501 503]));
mB = changeObjective(H7modelBioDemandnadh,H7rxnNamesNADH{1});[grRatioH7nadh grRateKOH7nadh grRateWTH7nadh] = singleGeneDeletion(mB);[fH7nadh,f_koH7nadh,vH7nadh,v_koH7nadh,status1H7nadh,lostxnsH7nadh,probtfgeneH7nadh] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);
mB = changeObjective(H7modelBioDemandnadh,H7rxnNamesNADH{2});[grRatioH7nadph grRateKOH7nadph grRateWTH7nadph] = singleGeneDeletion(mB);[fH7nadph,f_koH7nadph,vH7nadph,v_koH7nadph,status1H7nadph,lostxnsH7nadph,probtfgeneH7nadph] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);
[H7modelBioDemandnadh H7rxnNamesNADH] = addDemandReaction(H7modelDemand,H7modelDemand.mets([501 503]));
mB = changeObjective(H7modelBioDemandnadh,H7rxnNamesNADH{1});[grRatioH7nadh grRateKOH7nadh grRateWTH7nadh] = singleGeneDeletion(mB);[fH7nadh,f_koH7nadh,vH7nadh,v_koH7nadh,status1H7nadh,lostxnsH7nadh,probtfgeneH7nadh] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);
mB = changeObjective(H7modelBioDemandnadh,H7rxnNamesNADH{2});[grRatioH7nadph grRateKOH7nadph grRateWTH7nadph] = singleGeneDeletion(mB);[fH7nadph,f_koH7nadph,vH7nadph,v_koH7nadph,status1H7nadph,lostxnsH7nadph,probtfgeneH7nadph] =  promv2(mB,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);