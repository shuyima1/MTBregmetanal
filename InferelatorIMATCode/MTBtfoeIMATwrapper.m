% initCobraToolbox
% changeCobraSolver('glpk','MILP');

load MTBimatinput.mat Beste7H9 gens TFOEfccontrolExp

[B7IMATmodel, B7IMATrxns, mtbTFOEexpressionData] = expressionmatrixIMAT(Beste7H9,gens,TFOEfccontrolExp,-1,'ABSOLUTE','IMAT',0.001);

save MTBtfoeIMATresults B7* mtbTFOEexpressionData

