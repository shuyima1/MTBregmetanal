initCobraToolbox

%Finds the reactions in Fang that carry flux greater than 10^-10
[FangReduced hasFlux v1 v2] = reduceModel(iNJ661m,10^-10);
FanghasFlux = hasFlux;
Fangv1 = v1; %Vmax
Fangv2 = v2; %Vmin
soln = optimizeCbModel(FangReduced);
NoFluxMets= setdiff(iNJ661m.mets,FangReduced.mets);
FangNoFluxMets = NoFluxMets;
FangNoFluxRxns = setdiff(iNJ661m.rxns,FangReduced.rxns);

Fangsoln = optimizeCbModel(iNJ661m);
Fangsoln

%Find the reactions that produce NAD and ATP
nadx = strmatch('nad',iNJ661m.mets)
iNJ661m.mets(nadx)
atpx = strmatch('atp',iNJ661m.mets)
iNJ661m.mets(atpx)

nadix = cell(size(nadx,1),1)
for i = 1:size(nadx,1);nadix = find(iNJ661m.S(nadx(i),:));end
nadix = cell(size(nadx,1),1)
for i = 1:size(nadx,1);nadix{i} = find(iNJ661m.S(nadx(i),:) > 0);end

atpix = cell(size(atpx,1),1)
for i = 1:size(atpx,1);atpix{i} = find(iNJ661m.S(atpx(i),:) > 0);end
printRxnFormula(iNJ661m,iNJ661m.rxns(find(iNJ661m.S(atpx(1),:))))
printRxnFormula(iNJ661m,iNJ661m.rxns(find(iNJ661m.S(atpx(2),:))))
FangvATPc = [sum(Fangv1(atpix{1})) sum(Fangv2(atpix{1}))]
printRxnFormula(iNJ661m,iNJ661m.rxns(find(iNJ661m.S(atpx(2),:))))
find(iNJ661m.S(atpx(2),:)) 
iNJ661m.rxns(ans)
% Looks like atp[e] only has the exchange reaction coming into the system. 
% It never gets converted into atp[c] so that it can be used elsewhere.
% Therefore, the exchange flux on EX_atp(e) is approx zero.
Fangv1(266)
Fangv2(266)

% Comparing the exchange situation against oxygen, which has a reaction
% that can reversibly convert the o2[e] to o2[c] (O2t)

strmatch('EX_o',iNJ661m.rxns)
iNJ661m.rxns(ans)
strmatch('o2',iNJ661m.mets)
iNJ661m.mets(ans)
find(iNJ661m.S(551,:))
printRxnFormula(iNJ661m,iNJ661m.rxns(find(iNJ661m.S(551,:))))
iNJ661m.lb(723)
iNJ661m.ub(723)
Fangv2(723)
Fangv1(723)

% This would add the corresponding transfer reaction for atp
FangATPtr = addReaction(iNJ661m,'ATPt',[iNJ661m.mets(atpx(2));iNJ661m.mets(atpx(1))],[-1 1],true,-1000,1000)
soln2 = optimizeCbModel(FangATPtr)
soln2.x(1050)

% Examine more carefully the sources of atp[c]
printRxnFormula(iNJ661m,iNJ661m.rxns(atpix{1}))
FangvATPc = [Fangv1(atpix{1}(1)) Fangv2(atpix{1}(1));Fangv1(atpix{1}(2)) Fangv2(atpix{1}(2))]

% Finds which reactions (1) produce and (2) consume  (1) atp[c] and (2) atp[e]
atpix = cell(size(atpx,1),2)
for i = 1:size(atpx,1);atpix{i,1} = find(iNJ661m.S(atpx(i),:) > 0);atpix{i,2} = find(iNJ661m.S(atpx(i),:) < 0);end
atpR = iNJ661m.rev(atpix{1,2});
atpRix = find(atpR)
atpRfva = [Fangv1(atpix{1,2}(atpRix))' Fangv2(atpix{1,2}(atpRix))'];
FangvATPcAll = [sum(Fangv1([atpix{1,1} atpix{1,2}])) sum(Fangv2([atpix{1,1} atpix{1,2}]))]

% There are several reactions that have unconstrained FVA flux results for
% atp (ie vmin = -1000, vmax = 1000)
atpRix(25)
printRxnFormula(iNJ661m,iNJ661m.rxns(atpix{1,2}(atpRix([3 4 25]))))
