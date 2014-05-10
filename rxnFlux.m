function [flux] = rxnFlux(model,reactionID,varargin)

if nargin == 4
    soln = varargin{1};
    v1 = varargin{2}{1};
    v2 = varargin{2}{2};
elseif nargin == 2
    soln = optimizeCbModel(model,[],'one');
    [v1 v2] = fluxVariability(model);
end

id = strmatch(reactionID,model.rxns,'exact');
flux = printRxnFormula(model,reactionID);
flux = [reactionID flux num2cell([soln.x(id) v1(id) v2(id) model.lb(id) model.ub(id)])];

