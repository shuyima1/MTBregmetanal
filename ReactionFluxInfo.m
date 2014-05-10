function [MetRxns, MetRxnFormula, MetRxnIX, MetRxnFluxes] = ReactionFluxInfo(model,Met,FluxMatrix)

[MetRxns,MetRxnFormula] = findRxnsFromMets(model,Met);
MetRxnIX = findRxnIDs(model,MetRxns);
MetRxnFluxes = FluxMatrix(MetRxnIX,:);
