%Check ATP demand effect on Fang Aerobic Model
[FangmodelDemand FangrxnNames] = addDemandReaction(iNJ661m,iNJ661mDemandRResults.iNJ661mMets);
FangmodelDemand = changeObjective(FangmodelDemand,FangrxnNames{21});
[grRatioFangATP grRateKOFangATP grRateWTFangATP hasEffectFangATP delRxnsFangATP fluxSolnFangATP] = singleGeneDeletion(H1dbio);
[fFangATP,f_koFangATP,vFangATP,v_koFangATP,status1FangATP,lostxnsFangATP,probtfgeneFangATP] =  promv2(H1dbio,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);

%Check ATP demand effect on Fang Hypoxia Model
[H1modelDemand H1rxnNames] = addDemandReaction(H1model,iNJ661mDemandRResults.iNJ661mMets);
H1dbio = changeObjective(H1modelDemand,H1rxnNames{21});
[grRatioH1atp grRateKOH1atp grRateWTH1atp hasEffectH1atp delRxnsH1atp fluxSolnH1atp] = singleGeneDeletion(H1dbio);
[fH1atp,f_koH1atp,vH1atp,v_koH1atp,status1H1atp,lostxnsH1atp,probtfgeneH1atp] =  promv2(H1dbio,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);
fH1atp(38) % dosR is TF #38


%% Beste Model Analysis
%Constructing the Beste Hypoxia model
BesteO2_01 = changeRxnBounds(BesteVitroObjFuncBiomassE,'R804',0.01,'u'); % R804 is the oxygen exchange reaction
[BesteH1 BH1rxns] = createTissueSpecificModelSM(BesteO2_01,H1exp);
[grRatioBh1 grRateKOBh1 grRateWTBh1 hasEffectBh1 delRxnsBh1 fluxSolnBh1] = singleGeneDeletion(BesteH1);
[fH1B,f_koH1B,vH1B,v_koH1B,status1H1B,lostxnsH1B,probtfgeneH1B] =  promv2(BesteH1,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);

[BH1modelDemand BH1rxnNames] = addDemandReaction(BesteH1,BesteVitroDemandRResults.BesteVitroObjFuncBiomassEbiomassMets);

%Check TAG demand effect
BH1dTAGbio = changeObjective(BH1modelDemand,BH1rxnNames{3});
[grRatioH1Btag grRateKOH1Btag grRateWTH1Btag hasEffectH1Btag delRxnsH1Btag fluxSolnH1Btag] = singleGeneDeletion(BH1dTAGbio);
[fH1Btag,f_koH1Btag,vH1Btag,v_koH1Btag,status1H1Btag,lostxnsH1Btag,probtfgeneH1Btag] =  promv2(BH1dTAGbio,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);

%Check ATP demand effect
BH1dATPbio = changeObjective(BH1modelDemand,BH1rxnNames{4});
[grRatioBH1atp grRateKOBH1atp grRateWTBH1atp hasEffectBH1atp delRxnsBH1atp fluxSolnBH1atp] = singleGeneDeletion(BH1dATPbio);
[fH1Batp,f_koH1Batp,vH1Batp,v_koH1Batp,status1H1Batp,lostxnsH1Batp,probtfgeneH1Batp] =  promv2(BH1dATPbio,TFOEexp,TargetGene(:,1),BIGtf_05_2013(:,1),BIGtf_05_2013(:,2),[],[],[],[],[],[],[],[],0);

