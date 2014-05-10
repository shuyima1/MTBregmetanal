function [InferelatorIMATmodel, InferelatorIMATrxns, InferelatorKOexpressionData] = expressionmatrixIMAT(model,Genes,InferelatorKOexp,thresh,threshtype,runIMATflag,epsilon)

InferelatorIMATmodel = cell(size(InferelatorKOexp,2),1);
InferelatorIMATrxns = cell(size(InferelatorKOexp,2),1);

for nTF = 1:size(InferelatorKOexp,2)
    
   InferelatorKOexpressionData(nTF,1).Locus = Genes;
   if isequal(upper(threshtype),'QUANTILE')
       InferelatorKOexpressionData(nTF,1).Data = [InferelatorKOexp(:,nTF) > quantile(InferelatorKOexp(:,nTF),thresh)];
   elseif isequal(upper(threshtype),'ABSOLUTE')
       InferelatorKOexpressionData(nTF,1).Data = [InferelatorKOexp(:,nTF) > thresh];
   end
   InferelatorKOexpressionData(nTF,1).Data(isnan(InferelatorKOexp(:,nTF))) = 1;
   
   if isequal(upper(runIMATflag),'IMAT')
       % Run iMAT on the resulting model
       IMATmodel = [];
       IMATrxns = [];
       try
            [IMATmodel, IMATrxns] = createTissueSpecificModelimatorphanSM(model,InferelatorKOexpressionData(nTF,1),1,1,[],'iMAT',epsilon);
       catch err
           InferelatorIMATmodel{nTF} = nan;
           InferelatorIMATrxns{nTF} = nan;
       end
       InferelatorIMATmodel{nTF} = IMATmodel;
       InferelatorIMATrxns{nTF} = IMATrxns;
   elseif isequal(upper(runIMATflag),'GIMME')
       [IMATmodel, IMATrxns] = createTissueSpecificModelimatorphanSM(model,InferelatorKOexpressionData(nTF,1),1,1,[],'GIMME');
       InferelatorIMATmodel{nTF} = IMATmodel;
       InferelatorIMATrxns{nTF} = IMATrxns;
   end
end