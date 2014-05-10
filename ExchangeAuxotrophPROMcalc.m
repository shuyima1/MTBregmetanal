aafluxes = {'R820';'R821';'R822';'R823';'R826';'R829';
    'R830';'R831';'R832';'R833';'R834';'R835';'R839';
    'R842';'R843';'R848';'R849';'R850';'R864';'R866';'R868';'R869';
    'R870';'R871';'R874';'R875';'R877';'R878';'R881';'R883';'R884';'R886';
    'R887''R891';'R892';'R893';'R894';'R897';'R898'};

SautonFluxes = {'R800';'R804';'R805';'R812';'R822';'R841';'R851';'R882';'R924'};
% % add asparagine
% BesteSautonsJ = changeRxnBounds(Beste7H9modJ,'R822',1,'u');
% % block glucose
% BesteSautonsJ = changeRxnBounds(BesteSautonsJ,'R863',0,'u');
% % block glutamine
% BesteSautonsJ = changeRxnBounds(BesteSautonsJ,'R830',0,'u');
% BesteSautonsJ = changeRxnBounds(BesteSautonsJ,'R864',0,'u');
% % block Biotin
% BesteSautonsJ = changeRxnBounds(BesteSautonsJ,'R925',0,'u');


AuxExchanges = setdiff(BesteExchange(:,1),[aafluxes;SautonFluxes]);

AuxModelsSautons = cell(size(AuxExchanges,1),1);
for i = 1:size(AuxExchanges,1);
    AuxModelsSautons{i,1} = changeRxnBounds(BesteSautonsJ,AuxExchanges{i,1},1,'u');

    [grRatio grWT grWT delRs delRs sgdflux] = singleGeneDeletion(AuxModelsSautons{i,1});
    AuxModelsSautons{i,2} = {grRatio; grWT; delRs; sgdflux};
    [f fko v vko stat lostxns probtfgene] = promv2(AuxModelsSautons{i,1},metRvFC2,RvGeneIDs(:,1),BIGtf_05_2013_Rv2621(:,1),BIGtf_05_2013_Rv2621(:,2),[],[],[],[],[],[],0,[],0);
    AuxModelsSautons{i,3} = {f; v; fko; vko; probtfgene; lostxns; stat};
end

save auxotrophModelsmodJS AuxModelsSautons
