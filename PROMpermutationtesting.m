function [fnullexpperm, fnullnetperm] = PROMpermutationtesting(model,expression,expressionid,regulator,regulated,litevidence,prob_prior,subsets,v11,v12,KAPPA,datathresh,probtfgene,sizeflag,numperm,saveflag)

if (~exist('numperm','var')) || (isempty(numperm))
    numperm = 500;
end

if (~exist('v11','var')) || (isempty(v11))
        [v11,v12] = fluxVariability(model); 
end


fnullexpperm = zeros(size(unique(regulator(:,1)),1),numperm);
for i = 1:numperm
    expperm = expression(randperm(size(expression,1)),:);
    fnullexpperm(:,i) = promv2(model,expperm,expressionid(:,1),regulator,regulated,litevidence,prob_prior,subsets,v11,v12,KAPPA,datathresh,probtfgene,sizeflag);
end

if saveflag
    save PROMpermutation fnullexpperm
end

fnullnetperm = zeros(size(unique(regulator),1),numperm);
for i = 1:numperm
    netperm = [regulator regulated(randperm(size(regulated,1)))];
    fnullnetperm(:,i) = promv2(model,expression,expressionid(:,1),netperm(:,1),netperm(:,2),prob_prior,subsets,v11,v12,KAPPA,datathresh,probtfgene,sizeflag);
end

if saveflag
    save PROMpermutation fnullnetperm -append
end

end

    