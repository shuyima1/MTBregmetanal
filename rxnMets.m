function [posrxns negrxns soln v1 v2] = rxnMets(model,met,varargin)

if nargin == 4
    soln = varargin{1};
    v1 = varargin{2}{1};
    v2 = varargin{2}{2};
elseif nargin == 2
    soln = optimizeCbModel(model,[],'one');
    [v1 v2] = fluxVariability(model);
end

x = strmatch(met,model.mets,'exact');

posx = find(model.S(x,:) > 0);
negx = find(model.S(x,:) < 0);

disp('negative reactions:')
negrxns = printRxnFormula(model,model.rxns(negx));
negrxns = [model.rxns(negx) negrxns num2cell([soln.x(negx) v1(negx) v2(negx)])];

disp('positive reactions:')
posrxns = printRxnFormula(model,model.rxns(posx));
posrxns = [model.rxns(posx) posrxns num2cell([soln.x(posx) v1(posx) v2(posx)])];