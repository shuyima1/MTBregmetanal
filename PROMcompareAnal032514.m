initCobraToolbox
load('Jamshidimodels.mat', 'jamshidiBiGG')
load('BesteMediaModels.mat', 'Beste7H9')
load('phoPkoflux.mat')
load('BesteSautons.mat', 'BesteSautonsYesSulfate')
load('mtbPROMinputs2.mat', 'BIGtf_05_2013_Rv2621')
load('mtbPROMinputs2.mat', 'RvGeneIDs')
load('mtbPROMinputs2.mat', 'metRvFC2')

% WT fluxes and growth rates under different conditions
sB7 = optimizeCbModel(Beste7H9);
sSau = optimizeCbModel(BesteSautonsYesSulfate);
WTflux_7H9_Sauton = [sB7.x sSau.x];

[vminB7 vmaxB7] = fluxVariability(Beste7H9);
[vminSau vmaxSau] = fluxVariability(BesteSautonsYesSulfate);
WTfva_7H9_Sauton = [vminB7 vminSau vmaxB7 vmaxSau];

% PROM outputs
[fB7 fkoB7 vB7 vkoB7 statB7 lostxnsB7 probtfgeneB7] = promv2(Beste7H9,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fSau fkoSau vSau vkoSau statSau lostxnsSau probtfgeneSau] = promv2(BesteSautonsYesSulfate,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7pos fkoB7pos vB7pos vkoB7pos statB7 lostxnsB7 probtfgeneB7pos] = promv2pos(Beste7H9,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
fB7pos
[fSaupos fkoSaupos vSaupos vkoSaupos statSau lostxnsSau probtfgeneSaupos] = promv2pos(BesteSautonsYesSulfate,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);

uTF = unique(BIGtf_05_2013_Rv2621(:,1));
uTF{29} = 'Rv2621c';

% Rv0757-specific analyses
Rv0757flux_7H9_Sauton = [vB7(12,1:size(Beste7H9.rxns,1));vSau(12,1:size(Beste7H9.rxns,1))]';

phoPtargens
printRxnFormula(Beste7H9,'R925')
strmatch('R925',Beste7H9.rxns)
Beste7H9.lb(ans)
Beste7H9.ub(ans)
Beste7H9.ub(850)
BesteSautonsYesSulfate.ub(850)
BesteSautonsYesSulfate.lb(850)

size(probtfgeneB7)

xPhoPtf = strmatch('Rv0757',BIGtf_05_2013_Rv2621(:,1));
size(xPhoPtf)
phoPinteractions = BIGtf_05_2013_Rv2621(xPhoPtf,:);
phoPinteractions(:,3) = num2cell(probtfgeneB7(xPhoPtf));
phoPinteractions(:,4) = num2cell(probtfgeneB7pos(xPhoPtf));

intersect(strmatch('Rv0757',BIGtf_05_2013_Rv2621(:,1)),strmatch('Rv2121c',BIGtf_05_2013_Rv2621(:,2)))
probtfgeneB7(4374)
probtfgeneB7pos(4374)
probtfgeneSau(4374)
probtfgeneSaupos(4374)

intersect(strmatch('Rv0757',BIGtf_05_2013_Rv2621(:,1)),strmatch('Rv2392',BIGtf_05_2013_Rv2621(:,2)))
probtfgeneB7pos(4348)
probtfgeneB7(4348)
probtfgeneSau(4348)
probtfgeneSaupos(4348)

% Compare with the CDC1551 Rv0757 transposon differential expression
% dataset from PATRIC
phoPtransposonDE = readtext('TranscriptomicsGenes_CDC1551_Rv0757transposonVsWT.txt','\t');

[ia ia ib] = intersect(phoPinteractions(:,2),phoPtransposonDE(:,3),'stable');
phoPinteractions(ia,5:6) = phoPtransposonDE(ib,7:8);

load('TFOEsmall052813data.mat')
[ia2 ia2 ib2] = intersect(phoPinteractions(:,2),gens(:,1),'stable');
phoPinteractions(ia2,7) = tfoeExp(ib2,68);
phoPinteractions(ia2,7) = num2cell(tfoeExp(ib2,68));
phoPinteractions(ia2,8) = num2cell(tfoePval(ib2,68));

% regression correlation of expression between the PATRIC transposon data
% and the TFOE data
regcor = [];

[rho pv] = corr(regcor(~isnan(regcor(:,1)),1),regcor(~isnan(regcor(:,1)),2),'type','Spearman')
[rho pv] = corr(regcor(~isnan(regcor(:,1)),1),regcor(~isnan(regcor(:,1)),2))

[rho2 pv2] = corr(regcor(~isnan(regcor(:,1)),1),regcor(~isnan(regcor(:,1)),3))
[rho2 pv2] = corr(regcor(~isnan(regcor(:,1)),1),regcor(~isnan(regcor(:,1)),3),'type','Spearman')

[rho2 pv2] = corr(regcor(~isnan(regcor(:,1)),1),regcor(~isnan(regcor(:,1)),4),'type','Spearman')
[rho2 pv2] = corr(regcor(~isnan(regcor(:,1)),1),regcor(~isnan(regcor(:,1)),4))

[rho3 pv3] = corr(regcor(~isnan(regcor(:,1)),2),regcor(~isnan(regcor(:,1)),4))
[rho3 pv3] = corr(regcor(~isnan(regcor(:,1)),2),regcor(~isnan(regcor(:,1)),4),'type','Spearman')

sum(regcor(regcor(:,1) == -1,3) < 1)
size(regcor(regcor(:,1) == -1,3))
sum(regcor(regcor(:,1) == -1,3) == 1)
mean(regcor(regcor(:,1) == -1 & regcor(:,3) < 1,5))
mean(regcor(regcor(:,1) == -1 & regcor(:,3) == 1,5))
median(regcor(regcor(:,1) == -1 & regcor(:,3) == 1,5))
median(regcor(regcor(:,1) == -1 & regcor(:,3) < 1,5))

figure;hist(regcor(regcor(:,1) == -1 & regcor(:,3) == 1,5),100)
figure;hist(regcor(regcor(:,1) == -1 & regcor(:,3) < 1,5),100)

sum(regcor(regcor(:,1) == -1 & regcor(:,3) < 1,5) < 0.05)
sum(regcor(regcor(:,1) == -1 & regcor(:,3) == 1,5) < 0.05)

sum(regcor(:,1) == -1 & regcor(:,4) < 0)


% Single Gene Deletions of Metabolic genes in the different models
[grRatio7H9 grKO7H9 grWT7H9 h delRxns7H9 fluxKO7H9] = singleGeneDeletion(Beste7H9);
[grRatioSau grKOSau grWTSau h delRxnsSau fluxKOSau] = singleGeneDeletion(BesteSautonsYesSulfate);

find(grRatio7H9 < 0.01 & grRatioSau > 0.99)
Beste7H9.genes(ans)

flux2121KOcompare = [fluxKO7H9(:,378) fluxKOSau(:,378)];

% Need to change the flux of exchange reaction R831 (Histidine exchange) to zero
BesteSautonsYesSulfate.ub(756)
Beste7H9.ub(756)
BesteSautonsYesSulfate.ub(756) = 0;

load('BesteMediaModels.mat', 'BesteExchange')

%ia3 represents the indices of the exchange reactions in the models
[ia3 ia3 ib3] = intersect(BesteExchange(:,1),Beste7H9.rxns);

Bxch=[BesteExchange(ia3,1) num2cell([Beste7H9.ub(ib3) BesteSautonsYesSulfate.ub(ib3)])];

BesteSautonsYesSulfate = changeRxnBounds(BesteSautonsYesSulfate,'R811',0,'u');
BesteSautonsYesSulfate = changeRxnBounds(BesteSautonsYesSulfate,'R812',1,'u');

[grRatioSau grKOSau grWTSau h delRxnsSau fluxKOSau] = singleGeneDeletion(BesteSautonsYesSulfate);
[grRatio7H9 grKO7H9 grWT7H9 h delRxns7H9 fluxKO7H9] = singleGeneDeletion(Beste7H9);

cd gemini_files_code_data_0/

Bxch=[BesteExchange(ia3,1) num2cell([Beste7H9.ub(ib3) BesteSautonsYesSulfate.ub(ib3)])];

[fSau fkoSau vSau vkoSau statSau lostxnsSau probtfgeneSau] = promv2(BesteSautonsYesSulfate,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);

fSau = fSau';
fB7 = fB7';

find(grRatio7H9 > 0.99 & grRatioSau < 0.01)
Beste7H9.genes(ans)
find(grRatio7H9 < 0.01 & grRatioSau > 0.99)
Beste7H9.genes(ans)

sSau = optimizeCbModel(BesteSautonsYesSulfate)
sB7 = optimizeCbModel(Beste7H9)

cd ..
load BesteMediaModels BesteVitro

Bxch=[BesteExchange(ia3,1) num2cell([Beste7H9.ub(ib3) BesteSautonsYesSulfate.ub(ib3) BesteVitro.ub(ib3)])];
find(BesteVitro.c)
optimizeCbModel(BesteVitro)

% Since 7H9 medium contains bovine serum albumin, may need to open up
% fluxes of amino acids
Beste7H9aa = changeRxnBounds(Beste7H9,'R820',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R821',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R822',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R823',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R826',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R829',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R830',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R831',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R832',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R833',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R834',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R835',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R839',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R842',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R843',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R848',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R849',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R850',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R864',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R866',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R868',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R869',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R870',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R871',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R874',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R875',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R877',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R878',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R881',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R883',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R884',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R886',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R891',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R887',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R892',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R893',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R894',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R897',1,'u');
Beste7H9aa = changeRxnBounds(Beste7H9aa,'R904',1,'u');

findRxnsFromMets(Beste7H9aa,'BSA[c]')

Bxch=[BesteExchange(ia3,1) num2cell([Beste7H9.ub(ib3) BesteSautonsYesSulfate.ub(ib3) BesteVitro.ub(ib3) Beste7H9aa.ub(ib3)])];
optimizeCbModel(Beste7H9aa)

cd gemini_files_code_data_0/;
[fB7aapos fkoB7aapos vB7aapos vkoB7aapos statB7aa lostxnsB7aa probtfgeneB7aapos] = promv2pos(Beste7H9aa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7aa fkoB7aa vB7aa vkoB7aa statB7aa lostxnsB7aa probtfgeneB7aa] = promv2(Beste7H9aa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
fB7aa = fB7aa';

load Griffin/GriffinData.mat

[grRatio7H9aa grKO7H9aa grWT7H9aa h delRxns7H9aa fluxKO7H9aa] = singleGeneDeletion(Beste7H9);

[cG iaG ibG] = intersect(Beste7H9aa.genes,GriffinData(:,1),'stable');
cG(:,2:3) = GriffinData(ibG,7:8);
cG(:,4) = num2cell(grRatio7H9aa(iaG));

b1 = cell2mat(cG(:,2)) < 0.1;
b2 = strcmp('essential',cG(:,3));
b1a = cell2mat(cG(:,2)) > 0.9;
b2a = strcmp('non-essential',cG(:,3));

sum(cell2mat(cG(:,4)) > 0.95 & (b1a | b2a))
sum(cell2mat(cG(:,4)) > 0.95)

sum(cell2mat(cG(:,4)) < 0.05 & (b1 | b2))
sum(cell2mat(cG(:,4)) < 0.05)

149/174
303/499

%% compare the accuracies of the 7H9 models to sassetti 
cd gemini_files_code_data_0/
[fB7 fkoB7 vB7 vkoB7 statB7 lostxnsB7 probtfgeneB7] = promv2(Beste7H9,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7aa fkoB7aa vB7aa vkoB7aa statB7 lostxnsB7 probtfgeneB7aa] = promv2(Beste7H9aa,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);

fB7 = fB7';
fB7aa = fB7aa';

sB7a = optimizeCbModel(Beste7H9aa);
sB7 = optimizeCbModel(Beste7H9);

B7growthcompare=[fB7/sB7.f fB7aa/sB7a.f];

[grRatioB7 grKOB7 grWTB7 h delRxnsB7 fluxKOB7] = singleGeneDeletion(Beste7H9);
[grRatioB7aa grKOB7aa grWTB7aa h delRxnsB7aa fluxKOB7] = singleGeneDeletion(Beste7H9aa);

[cB7 iaB7 ibB7] = intersect(Beste7H9.genes,GriffinData(:,1),'stable');
cB7(:,2:3) = GriffinData(ibB7,7:8);
cB7(:,4) = num2cell(grRatioB7(iaB7));
cB7(:,5) = num2cell(grRatioB7aa(iaB7));
b2B7 = strcmp('essential',cB7(:,3));
b2B7a = strcmp('non-essential',cB7(:,3));

sum(cell2mat(cB7(:,4)) > 0.95 & (b2B7a))/sum(cell2mat(cB7(:,4)) > 0.95)
sum(cell2mat(cB7(:,4)) < 0.05 & (b2B7))/sum(cell2mat(cB7(:,4)) < 0.05)
sum(cell2mat(cB7(:,5)) > 0.95 & (b2B7a))/sum(cell2mat(cB7(:,5)) > 0.95)
sum(cell2mat(cB7(:,5)) < 0.05 & (b2B7))/sum(cell2mat(cB7(:,5)) < 0.05)

% looks like the Beste7H9 model without the amino acid fluxes does better
%%

[cGT iaGT ibGT] = intersect(uTF,GriffinData(:,1),'stable');
cGT(:,2:3) = GriffinData(ibGT,7:8);
cGT(:,4) = num2cell(fB7(iaGT)/sB7aa.f);
b1T = cell2mat(cGT(:,2)) < 0.1;
b2T = strcmp('essential',cGT(:,3));
b1Ta = cell2mat(cGT(:,2)) > 0.9;
b2Ta = strcmp('non-essential',cGT(:,3));
cGT(:,5) = num2cell(ones(size(cGT,1),1));
for i = 1:size(cGT,1);if b1T(i) | b2T(i); cGT{i,5} = 0;end;if b1Ta(i) | b2Ta(i); cGT{i,5} = 2;end;end
[r p] = corr(cell2mat(cGT(:,4)),cell2mat(cGT(:,5)))

[fB7 fkoB7 vB7 vkoB7 statB7 lostxnsB7 probtfgeneB7] = promv2(Beste7H9,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],-1,[],0);
fB7 = fB7';

cGT(:,4) = num2cell(fB7(iaGT)/sB7a.f);
[r p] = corr(cell2mat(cGT(:,4)),cell2mat(cGT(:,5)))



%% The Griffin media conditions are close to Sauton's but with addition of
% ethanol flux. Single metabolic Gene Deletion results improve with this
% flux change
BesteGriffin = changeRxnBounds(BesteSautonsYesSulfate,'R858',1,'u');
sBG = optimizeCbModel(BesteGriffin)

[grRatioG grKOG grWTG h delRxnsG fluxKOG] = singleGeneDeletion(BesteGriffin);
[cG2 iaG2 ibG2] = intersect(BesteGriffin.genes,GriffinData(:,1),'stable');
cG2(:,2:3) = GriffinData(ibG2,7:8);
cG2(:,4) = num2cell(grRatioG(iaG2));
b1G = cell2mat(cG2(:,2)) < 0.1;
b2G = strcmp('essential',cG2(:,3));
b1Ga = cell2mat(cG2(:,2)) > 0.9;
b2Ga = strcmp('non-essential',cG2(:,3));

sum(cell2mat(cG2(:,4)) > 0.95 & (b1Ga | b2Ga))
sum(cell2mat(cG2(:,4)) > 0.95)
sum(cell2mat(cG2(:,4)) < 0.05 & (b1G | b2G))
sum(cell2mat(cG2(:,4)) < 0.05)
195/224
294/456
%%

% run PROM on the Griffin media model
[fG fkoG vG vkoG statG lostxnsG probtfgeneG] = promv2(BesteGriffin,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
fG = fG';

% cG2T has the PROM results for the 50TF network with fold change
% experssion and the Griffin media metabolic model
[cG2T iaG2T ibG2T] = intersect(uTF,GriffinData(:,1),'stable');
cG2T(:,2:3) = GriffinData(ibG2T,7:8);
cG2T(:,4) = num2cell(fG(iaG2T)/sBG.f);
b1TG = cell2mat(cG2T(:,2)) < 0.1;
b2TG = strcmp('essential',cG2T(:,3));
b1TGa = cell2mat(cG2T(:,2)) > 0.9;
b2TGa = strcmp('non-essential',cG2T(:,3));
cG2T(:,5) = num2cell(ones(size(cG2T,1),1));
for i = 1:size(cG2T,1);if b1TG(i); cG2T{i,5} = 0;end;if b1TGa(i); cG2T{i,5} = 2;end;end

[r p] = corr(cell2mat(cG2T(:,4)),cell2mat(cG2T(:,5)),'type','Spearman')
[r p] = corr(cell2mat(cG2T(:,4)),cell2mat(cG2T(:,5)))
% PROM results don't have significant correlation with growth

ranksum(cell2mat(cG2T(cell2mat(cG2T(:,5)) == 0,4)),cell2mat(cG2T(cell2mat(cG2T(:,5)) == 2,4)))
mean(cell2mat(cG2T(cell2mat(cG2T(:,5)) == 0,4)))
mean(cell2mat(cG2T(cell2mat(cG2T(:,5)) == 2,4)))

[r p] = corr(cell2mat(cG2T(:,4)),cell2mat(cG2T(:,2)))

sum(cell2mat(cG2T(:,4)) < 0.99 & cell2mat(cG2T(:,2)) < 0.1)
sum(cell2mat(cG2T(:,4)) < 0.99)
sum(cell2mat(cG2T(:,4)) > 0.99 & cell2mat(cG2T(:,2)) > 0.99)
sum(cell2mat(cG2T(:,4)) > 0.99 & cell2mat(cG2T(:,2)) > 0.9)

sum(cell2mat(cG2T(:,4)) > 0.99 & cell2mat(cG2T(:,2)) > 0.1)
sum(cell2mat(cG2T(:,4)) < 0.5 & cell2mat(cG2T(:,2)) < 0.1)
sum(cell2mat(cG2T(:,4)) < 0.5 & cell2mat(cG2T(:,2)) < 0.25)
sum(cell2mat(cG2T(:,4)) < 0.5 & cell2mat(cG2T(:,2)) < 0.9)
sum(cell2mat(cG2T(:,4)) < 0.5)
sum(cell2mat(cG2T(:,4)) > 0.5 & cell2mat(cG2T(:,2)) > 0.9)
sum(cell2mat(cG2T(:,4)) > 0.5)

%% comparing the gene binarization threshold (from 0 to -1 log_2 fold change)
[fG fkoG vG vkoG statG lostxnsG probtfgeneG] = promv2(BesteGriffin,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],-1,[],0);
fG = fG';
cG2T(:,4) = num2cell(fG(iaG2T)/sBG.f);
[r p] = corr(cell2mat(cG2T(:,4)),cell2mat(cG2T(:,2)))
ranksum(cell2mat(cG2T(cell2mat(cG2T(:,5)) == 0,4)),cell2mat(cG2T(cell2mat(cG2T(:,5)) == 2,4)))
[r p] = corr(cell2mat(cG2T(:,4)),cell2mat(cG2T(:,5)))
% does not result in any change

%% try PROM on the current binding network. chipint is the ChIPseq binding network (03/18/14)
chipint = {};
chipint(1,:) = [];
chipint(strcmp('Rv2610c',chipint(:,1)),1) = {'Rv2621c'};
uTF2 = unique(chipint(:,1));

sBG = optimizeCbModel(BesteGriffin);

load tfoeindivarraysRvdata
RvexpFCmed = bsxfun(@minus,RvExp(:,16:end),RvExp(:,1));
[fG2 fkoG2 vG2 vkoG2 statG2 lostxnsG2 probtfgeneG2] = promv2(BesteGriffin,RvexpFCmed,Rvgens(:,1),chipint(:,1),chipint(:,2),[],[],[],[],[],[],0,[],0);
fG2 = fG2';

% cGbigT contains the PROM results for the large ChIPseq binding, all of
% the expression data (fold change) and the Griffin media metabolic model 
[cGbigT iaGbigT ibGbigT] = intersect(uTF2,GriffinData(:,1),'stable');
cGbigT(:,2:3) = GriffinData(ibGbigT,7:8);
sBG = optimizeCbModel(BesteGriffin);
cGbigT(:,4) = num2cell(fG2(iaGbigT)/sBG.f);
cGbigT{110,2} = nan;
cGbigT{101,2} = nan;
b1TbigG = cell2mat(cGbigT(:,2)) < 0.1;
b2TbigG = strcmp('essential',cGbigT(:,3));
b1TbigGa = cell2mat(cGbigT(:,2)) > 0.9;
b2TbigGa = strcmp('non-essential',cGbigT(:,3));
cGbigT(:,5) = num2cell(ones(size(cGbigT,1),1));
for i = 1:size(cGbigT,1);if b1TbigG(i); cGbigT{i,5} = 0;end;if b1TbigGa(i); cGbigT{i,5} = 2;end;end
[r p] = corr(cell2mat(cGbigT(:,4)),cell2mat(cGbigT(:,5)),'type','Spearman')
[r p] = corr(cell2mat(cGbigT(:,4)),cell2mat(cGbigT(:,5)))


%% try PROM on a tf network based only on significant TFOE expression
% changes (sigtfoe comes from the tfoesmall expression network
sigtfoe = abs(tfoeExp) > 1 & tfoePval < 0.01;
[sigrow sigcol] = find(sigtfoe);
sigrow = sigrow(sigcol > 11);
sigcol = sigcol(sigcol > 11);
tfoesigints = [samples(sigcol) gens(sigrow)];
uTF3 = unique(tfoesigints(:,1));

[fGtfoe fkoGtfoe vGtfoe vkoGtfoe statGtfoe lostxnsGtfoe probtfgeneGtfoe] = promv2(BesteGriffin,tfoeExp,gens(:,1),tfoesigints(:,1),tfoesigints(:,2),[],[],[],[],[],[],0,[],0);

% cGtfoeT contains the PROM results for a TFOE-based network,
% the median fold change expression data for all genes and the Griffin media metabolic model 
[cGtfoeT iaGtfoeT ibGtfoeT] = intersect(uTF3,GriffinData(:,1),'stable');
cGtfoeT(:,2:3) = GriffinData(ibGtfoeT,7:8);
cGtfoeT{128,2} = nan;
cGtfoeT{140,2} = nan;
cGtfoeT(:,4) = num2cell(fGtfoe(iaGtfoeT)/sBG.f);
b1TtfoeG = cell2mat(cGtfoeT(:,2)) < 0.1;
b2TtfoeG = strcmp('essential',cGtfoeT(:,3));
b1TtfoeGa = cell2mat(cGtfoeT(:,2)) > 0.9;
b2TtfoeGa = strcmp('non-essential',cGtfoeT(:,3));
cGtfoeT(:,5) = num2cell(ones(size(cGtfoeT,1),1));
for i = 1:size(cGtfoeT,1);if b1TtfoeG(i); cGtfoeT{i,5} = 0;end;if b1TtfoeGa(i); cGtfoeT{i,5} = 2;end;end
[r p] = corr(cell2mat(cGtfoeT(:,4)),cell2mat(cGtfoeT(:,5)),'type','Spearman')
[r p] = corr(cell2mat(cGtfoeT(:,4)),cell2mat(cGtfoeT(:,5)))

strmatch('Rv3583c',samples)
strmatch('Rv3583c',uTF3)
find(sigrow == 205)

%% try PROM with using an absolute value of gene expression rather than fold
% change
[fGchipabs fkoGchipabs vGchipabs vkoGchipabs statG2 lostxnsG2 probtfgeneGchipabs] = promv2(BesteGriffin,RvExp,Rvgens(:,1),chipint(:,1),chipint(:,2),[],[],[],[],[],[],quantile(RvExp(:),0.15),[],0);
cGbigTabs = cGbigT;
cGbigT(:,4) = num2cell(fGchipabs(iaGbigT)/sBG.f);
cGbigTabs(:,4) = num2cell(fGchipabs(iaGbigT)/sBG.f);
cGbigT(:,4) = num2cell(fG2(iaGbigT)/sBG.f);
[r p] = corr(cell2mat(cGbigT(:,4)),cell2mat(cGbigT(:,5)),'type','Spearman')
[r p] = corr(cell2mat(cGbigT(:,4)),cell2mat(cGbigT(:,5)))

%not as effective 

[c1 ia1 ib1] = intersect(cG2T(:,1),uTF2);
[c2 ia2 ib2] = intersect(cG2T(:,1),uTF3);

cG2T(ia1,6) = num2cell(fG2(ib1)/sBG.f);
cG2T(ia2,7) = num2cell(fGtfoe(ib2)'/sBG.f);

cG2T(ia1,8) = num2cell(fG2(ib1)/sBG.f);
cG2T(ia1,9) = num2cell(fGchipabs(ib1)'/sBG.f);

[fGchipabs33 fkoGchipabs33 vGchipabs33 vkoGchipabs33 statG2 lostxnsG2 probtfgeneGchipabs33] = promv2(BesteGriffin,RvExp,Rvgens(:,1),chipint(:,1),chipint(:,2),[],[],[],[],[],[],quantile(RvExp(:),0.33),[],0);
cG2T(ia1,10) = num2cell(fGchipabs33(ib1)'/sBG.f);
isequal(cG2T(:,9),cG2T(:,10))
isequal(probtfgeneGchipabs,probtfgeneGchipabs33)
[fGchipabs33 fkoGchipabs33 vGchipabs33 vkoGchipabs33 statG2 lostxnsG2 probtfgeneGchipabs33] = promv2(BesteGriffin,RvExp,Rvgens(:,1),chipint(:,1),chipint(:,2),[],[],[],[],[],[],quantile(RvExp(:),0.5),[],0);
isequal(probtfgeneGchipabs,probtfgeneGchipabs33)

%% comparing the PROM performances on the 50 TF set
[r p] = corr(cell2mat(cG2T(:,4)),cell2mat(cG2T(:,2)))
[r p] = corr(cell2mat(cG2T(:,6)),cell2mat(cG2T(:,2)))
[r p] = corr(cell2mat(cG2T(:,7)),cell2mat(cG2T(:,2)))
[r p] = corr(cell2mat(cG2T(:,8)),cell2mat(cG2T(:,2)))
[r p] = corr(cell2mat(cG2T(:,9)),cell2mat(cG2T(:,2)))
[r p] = corr(cell2mat(cG2T(:,10)),cell2mat(cG2T(:,2)))

%% try PROM with the small network flag
[fG1 fkoG1 vG1 vkoG1 statG1 lostxnsG1 probtfgeneG1] = promv2(BesteGriffin,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],1);
cG2T(:,10) = num2cell(fG1'/sBG.f);
isequal(cG2T(:,10),cG2T(:,4))
[r p] = corr(cell2mat(cG2T(:,10)),cell2mat(cG2T(:,2)))
[r p] = corr(cell2mat(cG2T(:,4)),cell2mat(cG2T(:,5)))
b1TG = cell2mat(cG2T(:,2)) < 0.15;
b2TG = strcmp('essential',cG2T(:,3));
b1TGa = cell2mat(cG2T(:,2)) > 0.85;
b2TGa = strcmp('non-essential',cG2T(:,3));
cG2T(:,5) = num2cell(ones(size(cG2T,1),1));
for i = 1:size(cG2T,1);if b1TG(i) | b2TG(i); cG2T{i,5} = 0;end;if b1TGa(i) | b2TGa(i); cG2T{i,5} = 2;end;end
[r p] = corr(cell2mat(cG2T(:,4)),cell2mat(cG2T(:,5)))
% not as good as with the large network flag

%% Test if there is a significant difference between the growth rates of essential, non-essential genes
ranksum(cell2mat(cG2T(cell2mat(cG2T(:,5)) == 0,4)),cell2mat(cG2T(cell2mat(cG2T(:,5)) > 0,4)))

ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.2,4)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.2,4)))
ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.15,4)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.15,4)))
ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.10,4)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.10,4)))

ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.15,6)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.15,6)))
ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.15,7)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.15,7)))

ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.15,8)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.15,8)))
ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.2,8)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.2,8)))
ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.1,8)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.1,8)))
ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.05,8)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.05,8)))

ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.1,9)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.1,9)))

ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.1,10)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.1,10)))
ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.15,10)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.15,10)))
ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.05,10)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.05,10)))

[r p] = corr(cell2mat(cGbigT(:,4)),cell2mat(cGbigT(:,2)),'type','Spearman')
[r p] = corr(cell2mat(cGbigT(:,4)),cell2mat(cGbigT(:,2)))

[r p] = corr(cell2mat(cGbigT(~isnan(cell2mat(cGbigT(:,2))),4)),cell2mat(cGbigT(~isnan(cell2mat(cGbigT(:,2))),2)),'type','Spearman')
[r p] = corr(cell2mat(cGbigT(~isnan(cell2mat(cGbigT(:,2))),4)),cell2mat(cGbigT(~isnan(cell2mat(cGbigT(:,2))),2)))

[r p] = corr(cell2mat(cGbigTabs(~isnan(cell2mat(cGbigTabs(:,2))),4)),cell2mat(cGbigTabs(~isnan(cell2mat(cGbigTabs(:,2))),2)),'type','Spearman')
[r p] = corr(cell2mat(cGbigTabs(~isnan(cell2mat(cGbigTabs(:,2))),4)),cell2mat(cGbigTabs(~isnan(cell2mat(cGbigTabs(:,2))),2)))

% try to compare the significance for a range of threshold values
blah = zeros(20,2);
for i = 1:20
    try
        blah(i,1) = ranksum(cell2mat(cGbigT(cell2mat(cGbigT(:,2)) < 0.05*i,4)),cell2mat(cGbigT(cell2mat(cGbigT(:,2)) > 0.05*i,4)));
    catch exception
        blah(i,1) = nan;
    end
end
for i = 1:20
    try
        blah(i,2) = ranksum(cell2mat(cGbigT(cell2mat(cGbigT(:,4)) < 0.05*i,2)),cell2mat(cGbigT(cell2mat(cGbigT(:,4)) > 0.05*i,2)));
    catch exception
        blah(i,2) = nan;
    end
end

blah2 = zeros(20,2);
for i = 1:20
    try
        blah2(i,2) = ranksum(cell2mat(cG2T(cell2mat(cG2T(:,4)) < 0.05*i,2)),cell2mat(cG2T(cell2mat(cG2T(:,4)) > 0.05*i,2)));
    catch exception
        blah2(i,2) = nan;
    end
end
for i = 1:20
    try
        blah2(i,1) = ranksum(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.05*i,4)),cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.05*i,4)));
    catch exception;
        blah2(i,1) = nan;
    end
end

%% Compare with the original PROM results
PROMoldresults = {};
[cp iap ibp] = intersect(PROMoldresults(:,1),GriffinData(:,1),'stable');
PROMoldresults(iap,3:4) = GriffinData(ibp,7:8);

blah3 = zeros(20,2);for i = 1:20;try;blah3(i,2) = ranksum(cell2mat(PROMoldresults(cell2mat(PROMoldresults(:,2))/0.052 < 0.05*i,3)),cell2mat(PROMoldresults(cell2mat(PROMoldresults(:,2))/0.052 > 0.05*i,3)));catch exception; blah3(i,2) = nan;end;end
for i = 1:20;try;blah3(i,1) = ranksum(cell2mat(PROMoldresults(cell2mat(PROMoldresults(:,3)) < 0.05*i,2)),cell2mat(PROMoldresults(cell2mat(PROMoldresults(:,3)) > 0.05*i,2)));catch exception; blah3(i,1) = nan;end;end

[r p] = corr(cell2mat(PROMoldresults(:,2))/0.052,cell2mat(PROMoldresults(:,3)))

cell2mat(PROMoldresults(:,5))
PROMoldresults(:,5) = num2cell(strcmp('essential',PROMoldresults(:,4)));
ranksum(cell2mat(PROMoldresults(cell2mat(PROMoldresults(:,5)),2)),cell2mat(PROMoldresults(~cell2mat(PROMoldresults(:,5)),2)))
%% Try to compare how well the PROM results stack up against just Sassetti

ranksum(cell2mat(cGbigT(strcmp('essential',cGbigT(:,3)),4)),cell2mat(cGbigT(~strcmp('essential',cGbigT(:,3)),4)))
ranksum(cell2mat(cG2T(strcmp('essential',cG2T(:,3)),4)),cell2mat(cG2T(~strcmp('essential',cG2T(:,4)),2)))

%%
cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.15,4))
cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.15,4))
mean(cell2mat(cG2T(cell2mat(cG2T(:,2)) > 0.15,4)))
mean(cell2mat(cG2T(cell2mat(cG2T(:,2)) < 0.15,4)))

ranksum(cell2mat(cGbigT(cell2mat(cGbigT(:,4)) < 0.99,2)),cell2mat(cGbigT(cell2mat(cGbigT(:,4)) > 0.99,2)))

[r p] = corr(cell2mat(cG2T(:,4)),cell2mat(cG2T(:,2)))
[r p] = corr(cell2mat(cG2T(:,4)),cell2mat(cG2T(:,2)),'type','Spearman')

load('BesteMediaModels.mat', 'Beste7H9')

Bxch=[BesteExchange(ia3,1) num2cell([Beste7H9.ub(ib3) BesteSautonsYesSulfate.ub(ib3) BesteVitro.ub(ib3) Beste7H9aa.ub(ib3)])];

