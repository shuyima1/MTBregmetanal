% Beste7H9Jaa2 = changeRxnBounds(Beste7H9modJ,'R820',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R821',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R822',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R823',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R826',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R829',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R830',1,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R831',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R832',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R833',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R834',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R835',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R839',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R842',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R843',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R848',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R849',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R850',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R864',1,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R866',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R868',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R869',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R870',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R871',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R874',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R875',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R877',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R878',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R881',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R883',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R884',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R886',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R891',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R887',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R892',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R893',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R894',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R897',10^-4,'u');
% Beste7H9Jaa2 = changeRxnBounds(Beste7H9Jaa2,'R898',10^-4,'u');

aaXR = cell(20,1);
aaXR{1} = 'R820'; % ala
aaXR{2} = 'R821'; % arg
aaXR{3} = 'R822'; % asn
aaXR{4} = {'R823';'R848';'R849';'R850'}; % asp
aaXR{5} = 'R826'; % cys
aaXR{6} = 'R829'; % gln
aaXR{7} = {'R830';'R864'}; % glu
aaXR{8} = 'R866'; %gly
aaXR{9} = {'R831';'R868';'R869'}; % his
aaXR{10} = {'R832';'R870';'R871'}; % ile
aaXR{11} = {'R833';'R874';'R875'}; % leu
aaXR{12} = {'R834';'R877';'R878'}; % lys
aaXR{13} = 'R835'; % met
aaXR{14} = 'R881'; % phe
aaXR{15} = {'R839';'R883';'R884'}; % pro
aaXR{16} = {'R886';'R887'}; % ser
aaXR{17} = {'R842';'R891';'R892'}; % thr
aaXR{18} = 'R893'; % trp
aaXR{19} = 'R894'; % tyr
aaXR{20} = {'R843';'R897';'R898'}; % val

aaModelsSautons = cell(20,1);
for i = 1:20;
    aaModelsSautons{i,1} = changeRxnBounds(BesteSautonsJ,aaXR{i,1},1,'u');
%     aaModels{i,1} = changeRxnBounds(Beste7H9modJ,aaXR{i,1},1,'u');
end

for i = 1:20
    %aaModels{i,2} = optimizeCbModel(aaModels{i,1});
    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(aaModelsSautons{i,1});
    aaModelsSautons{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(aaModelsSautons{i,1},metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
    aaModelsSautons{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

save aaModelsmodJS aaModelsSautons