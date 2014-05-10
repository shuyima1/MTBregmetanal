[grRatio grKO grWT] = singleGeneDeletion(Beste7H9);
[grRatioJ grKOJ grWTJ] = singleGeneDeletion(Beste7H9modJ);
[c1 ia1 ib1] = intersect(Beste7H9.genes,GriffinData(:,1));
[c2 ia2 ib2] = intersect(Beste7H9modJ.genes,GriffinData(:,1));

c1(:,2:3) = GriffinData(ib1,7:8);
c1(:,4) = num2cell(grRatio(ia1));
c2(:,2:3) = GriffinData(ib2,7:8);
c2(:,4) = num2cell(grRatioJ(ia2));

sum(cell2mat(c1(strcmp('essential',c1(:,3))|strcmp('growth-defect',c1(:,3)),4)) < 0.85)/sum(strcmp('essential',c1(:,3))|strcmp('growth-defect',c1(:,3)))
%0.6160
sum(cell2mat(c1(strcmp('non-essential',c1(:,3)),4)) > 0.85)/sum(strcmp('non-essential',c1(:,3)))
%0.8338
% average of the two is: 0.7249
B7H9sgdperf = c1;

sum(cell2mat(c2(strcmp('essential',c2(:,3))|strcmp('growth-defect',c2(:,3)),4)) < 0.85)/sum(strcmp('essential',c2(:,3))|strcmp('growth-defect',c2(:,3)))
%0.5824
sum(cell2mat(c2(strcmp('non-essential',c2(:,3)),4)) > 0.85)/sum(strcmp('non-essential',c2(:,3)))
%0.8412
% average of the two is: 0.7118

[ce iae ibe] = intersect(BesteExchange(:,1),Beste7H9.rxns,'stable');
[ce2 iae2 ibe2] = intersect(BesteExchange(:,1),Beste7H9modJ.rxns,'stable');

% add some amino acid flux from source: BSA
Beste7H9Jaa2 = changeRxnBounds(Beste7H9modJ,'R820',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R821',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R822',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R823',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R826',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R829',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R830',1,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R831',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R832',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R833',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R834',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R835',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R839',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R842',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R843',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R848',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R849',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R850',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R864',1,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R866',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R868',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R869',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R870',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R871',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R874',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R875',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R877',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R878',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R881',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R883',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R884',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R886',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R891',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R887',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R892',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R893',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R894',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R897',10^-4,'u');
Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R898',10^-4,'u');

% add asparagine
BesteGriffinJ = changeRxnBounds(Beste7H9modJ,'R822',1,'u');
% add ethanol
BesteGriffinJ = changeRxnBounds(BesteGriffinJ,'R858',1,'u');
% block glucose
BesteGriffinJ = changeRxnBounds(BesteGriffinJ,'R863',0,'u');
% block glutamine
BesteGriffinJ = changeRxnBounds(BesteGriffinJ,'R830',0,'u');
BesteGriffinJ = changeRxnBounds(BesteGriffinJ,'R864',0,'u');
% block Biotin
BesteGriffinJ = changeRxnBounds(BesteGriffinJ,'R925',0,'u');

% add asparagine
BesteSautonsJ = changeRxnBounds(Beste7H9modJ,'R822',1,'u');
% block glucose
BesteSautonsJ = changeRxnBounds(BesteSautonsJ,'R863',0,'u');
% block glutamine
BesteSautonsJ = changeRxnBounds(BesteSautonsJ,'R830',0,'u');
BesteSautonsJ = changeRxnBounds(BesteSautonsJ,'R864',0,'u');
% block Biotin
BesteSautonsJ = changeRxnBounds(BesteSautonsJ,'R925',0,'u');

[grRatioJaa2 grKOJaa2 grWTJaa2] = singleGeneDeletion(Beste7H9Jaa2);
[grRatioJG grKOJG grWTJG] = singleGeneDeletion(BesteGriffinJ);
[grRatioJS grKOJS grWTJS] = singleGeneDeletion(BesteSautonsJ);

[c3 ia3 ib3] = intersect(Beste7H9Jaa2.genes,GriffinData(:,1));
c3(:,2:3) = GriffinData(ib3,7:8);
c3(:,4) = num2cell(grRatioJaa2(ia3));
c3(:,5) = num2cell(grRatioJG(ia3));
c3(:,6) = num2cell(grRatioJS(ia3));

sum(cell2mat(c3(strcmp('essential',c3(:,3))|strcmp('growth-defect',c3(:,3)),4)) < 0.85)/sum(strcmp('essential',c3(:,3))|strcmp('growth-defect',c3(:,3)))
% 0.5824
sum(cell2mat(c3(strcmp('non-essential',c3(:,3)),4)) > 0.85)/sum(strcmp('non-essential',c3(:,3)))
% 0.8412
% average of the two is: 0.7184

sum(cell2mat(c3(cell2mat(c3(:,2)) < 0.15,5)) < 0.85)/sum(cell2mat(c3(:,2)) < 0.15)
% 0.5423
sum(cell2mat(c3(cell2mat(c3(:,2)) > 0.15,5)) > 0.85)/sum(cell2mat(c3(:,2)) > 0.15)
% 0.8945
% average of the two is: 0.7369

load('mtbPROMinputs2.mat', 'BIGtf_05_2013_Rv2621')
load('mtbPROMinputs2.mat', 'metRvFC2')
load('mtbPROMinputs2.mat', 'RvGeneIDs')

[fB7Jaa2 fkoB7Jaa2 vB7Jaa2 vkoB7Jaa2 statB7 lostxnsB7 probtfgeneB7Jaa2] = promv2(Beste7H9Jaa2,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7JG fkoB7JG vB7JG vkoB7JG statB7 lostxnsB7 probtfgeneB7JG] = promv2(BesteGriffinJ,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7JS fkoB7JS vB7JS vkoB7JS statB7 lostxnsB7 probtfgeneB7JS] = promv2(BesteSautonsJ,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
uTF = unique(BIGtf_05_2013_Rv2621(:,1));
[cTF iaTF ibTF] = intersect(uTF,GriffinData(:,1));
cTF(:,2:3) = GriffinData(ibTF,7:8);
cTF(:,4) = num2cell(fB7Jaa2/grWTJaa2);
cTF(:,5) = num2cell(fB7JG/grWTJG);
cTF(:,6) = num2cell(fB7JS/grWTJS);

sum(cell2mat(cTF(strcmp('essential',cTF(:,3))|strcmp('growth-defect',cTF(:,3)),4)) < 0.85)/sum(strcmp('essential',cTF(:,3))|strcmp('growth-defect',cTF(:,3)))
% 1
sum(cell2mat(cTF(strcmp('non-essential',cTF(:,3)),4)) > 0.85)/sum(strcmp('non-essential',cTF(:,3)))
% 0.6410
% average of the two is: 0.8205

sum(cell2mat(cTF(cell2mat(cTF(:,2)) < 0.15,5)) < 0.85)/sum(cell2mat(cTF(:,2)) < 0.15)
% 0.7778
sum(cell2mat(cTF(cell2mat(cTF(:,2)) > 0.15,5)) > 0.85)/sum(cell2mat(cTF(:,2)) > 0.15)
% 0.6829
% average of the two is: 0.7304

%% acetate model
%acetate
Beste7H9Jac = changeRxnBounds(Beste7H9Jaa2,'R808',1,'u');
Beste7H9Jac = changeRxnBounds(Beste7H9Jac,'R844',1,'u');
%glycerol
Beste7H9Jac = changeRxnBounds(Beste7H9Jac,'R812',0,'u');
%glucose
Beste7H9Jac = changeRxnBounds(Beste7H9Jac,'R863',0,'u');

%% glucose model
Beste7H9Jglc = changeRxnBounds(Beste7H9Jaa2,'R812',0,'u');

%% glycerol model
Beste7H9Jglyc = changeRxnBounds(Beste7H9Jaa2,'R863',0,'u');

Bxt = [BesteExchange(iae,1:2) num2cell([Beste7H9modJ.ub(ibe2) BesteGriffinJ.ub(ibe2) BesteSautonsJ.ub(ibe2) Beste7H9Jaa2.ub(ibe2) Beste7H9Jac.ub(ibe2) Beste7H9Jglc.ub(ibe2) Beste7H9Jglyc.ub(ibe2)])];
[fB7Jac fkoB7Jac vB7Jac vkoB7Jac statB7 lostxnsB7 probtfgeneB7Jaa2] = promv2(Beste7H9Jac,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7Jglc fkoB7Jglc vB7Jglc vkoB7Jglc statB7 lostxnsB7 probtfgeneB7Jaa2] = promv2(Beste7H9Jglc,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7Jglyc fkoB7Jglyc vB7Jglyc vkoB7Jglyc statB7 lostxnsB7 probtfgeneB7Jaa2] = promv2(Beste7H9Jglyc,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);

[grRatioJac grKOJac grWTJac] = singleGeneDeletion(Beste7H9Jac);
[grRatioJglc grKOJglc grWTJglc] = singleGeneDeletion(Beste7H9Jglc);
[grRatioJglyc grKOJglyc grWTJglyc] = singleGeneDeletion(Beste7H9Jglyc);

cTF(:,7) = num2cell(fB7Jac/grWTJac);
cTF(:,8) = num2cell(fB7Jglc/grWTJglc);
cTF(:,9) = num2cell(fB7Jglyc/grWTJglyc);

%% acetate model
%acetate
Beste7H9JSac = changeRxnBounds(BesteSautonsJ,'R808',1,'u');
Beste7H9JSac = changeRxnBounds(Beste7H9JSac,'R844',1,'u');
%glycerol
Beste7H9JSac = changeRxnBounds(Beste7H9JSac,'R812',0,'u');
%glucose
Beste7H9JSac = changeRxnBounds(Beste7H9JSac,'R863',0,'u');

%% glucose model
Beste7H9Jglc = changeRxnBounds(BesteSautonsJ,'R812',0,'u');

%%
Bxt = [BesteExchange(iae,1:2) num2cell([Beste7H9modJ.ub(ibe2) BesteGriffinJ.ub(ibe2) BesteSautonsJ.ub(ibe2) Beste7H9Jaa2.ub(ibe2) Beste7H9Jac.ub(ibe2) Beste7H9Jglc.ub(ibe2) Beste7H9Jglyc.ub(ibe2)])];
[fB7Jac fkoB7Jac vB7Jac vkoB7Jac statB7 lostxnsB7 probtfgeneB7Jaa2] = promv2(Beste7H9Jac,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7Jglc fkoB7Jglc vB7Jglc vkoB7Jglc statB7 lostxnsB7 probtfgeneB7Jaa2] = promv2(Beste7H9Jglc,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
[fB7Jglyc fkoB7Jglyc vB7Jglyc vkoB7Jglyc statB7 lostxnsB7 probtfgeneB7Jaa2] = promv2(Beste7H9Jglyc,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);

[grRatioJac grKOJac grWTJac] = singleGeneDeletion(Beste7H9Jac);
[grRatioJglc grKOJglc grWTJglc] = singleGeneDeletion(Beste7H9Jglc);
[grRatioJglyc grKOJglyc grWTJglyc] = singleGeneDeletion(Beste7H9Jglyc);

cTF(:,7) = num2cell(fB7Jac/grWTJac);
cTF(:,8) = num2cell(fB7Jglc/grWTJglc);
cTF(:,9) = num2cell(fB7Jglyc/grWTJglyc);
cTF(:,10) = num2cell(fB7JSac/grWTJSac);
cTF(:,11) = num2cell(fB7JSglc/grWTJSglc);


aaModelsSimulations
fdoubleJaa2 = PROMdoubleKO(Beste7H9Jaa2,metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],0,[],0);save Beste7H9Jaa2DoubleKOresults fdoubleJaa2

aafluxes = {'R820';'R821';'R822';'R823';'R826';'R829';
'R830';'R831';'R832';'R833';'R834';'R835';'R839';
'R842';'R843';'R848';'R849';'R850';'R864';'R866';'R868';'R869';
'R870';'R871';'R874';'R875';'R877';'R878';'R881';'R883';'R884';'R886';
'R887''R891';'R892';'R893';'R894';'R897';'R898'};
SautonFluxes = {'R800';'R804';'R805';'R812';'R822';'R841';'R851';'R882';'R924'};
AuxExchanges = setdiff(Bxt(:,1),[aafluxes;SautonFluxes]);
size(AuxExchanges)



AuxModelsSautons = cell(size(AuxExchanges,1),1);
for i = 1:size(AuxExchanges,1);
    AuxModelsSautons{i,1} = changeRxnBounds(BesteSautonsJ,AuxExchanges{i,1},1,'u');
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(AuxModelsSautons{i,1});
    AuxModelsSautons{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(AuxModelsSautons{i,1},metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
    AuxModelsSautons{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end
save auxotrophModelsmodJS AuxModelsSautons

AuxJSaugrRatio = zeros(50,size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJSaugrRatio(:,i) = AuxModelsSautons{i,3}{1}/AuxModelsSautons{i,2}{2};end
AuxJSaugrKO = zeros(50,size(AuxExchanges,1));for i = 1:size(AuxExchanges,1);AuxJSaugrKO(:,i) = AuxModelsSautons{i,3}{1};end
AuxExchanges(:,2) = printRxnFormula(BesteSautonsJ,AuxExchanges(:,1));

ExtractTFinfoJsauton
TFrxnLongJS = cell(50,1);
for i = 1:50;eval(['flag = ~isempty(' uTF{i} 'affectedix);']);if flag;eval(['TFrxnLongJS{i} = ' uTF{i} 'RxnLongJS;']);end;end
for i = 1:50; if ~isempty(TFrxnLongJS{i});[ix ix] = unique(TFrxnLongJS{i}(:,1));TFrxnLongJS{i} = TFrxnLongJS{i}(ix,:);end;end

for i = 1:size(uTF,1)
    if ~isempty(TFrxnLongJS{i,1})
        [iax,iax,ibx] = intersect(TFrxnLongJS{i,1}(:,1),BesteSautonsJ.rxns,'stable');
        flag = 0;
        for j = 1:size(iax,1)
            if vB7JS(i,ibx(j)) >= 0
                if flag == 0
                    TFrxnLongJS{i,2} = changeRxnBounds(BesteSautonsJ,BesteSautonsJ.rxns(ibx(j)),vB7JS(i,ibx(j)),'u');
                    flag = 1;
                elseif flag == 1
                    TFrxnLongJS{i,2} = changeRxnBounds(TFrxnLongJS{i,2},BesteSautonsJ.rxns(ibx(j)),vB7JS(i,ibx(j)),'u');
                end
            elseif vB7JS(i,ibx(j)) < 0
                if flag == 0
                    TFrxnLongJS{i,2} = changeRxnBounds(BesteSautonsJ,BesteSautonsJ.rxns(ibx(j)),vB7JS(i,ibx(j)),'l');
                    flag = 1;
                elseif flag == 1
                    TFrxnLongJS{i,2} = changeRxnBounds(TFrxnLongJS{i,2},BesteSautonsJ.rxns(ibx(j)),vB7JS(i,ibx(j)),'l');
                end
            end
        end
    end
end


