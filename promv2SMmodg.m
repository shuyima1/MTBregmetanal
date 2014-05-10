function [f,f_ko,v,v_ko,status1,lostxns,probtfgene] =  promv2SMmodg(model,expression,expressionid,regulator,targets,litevidence,prob_prior,subsets,v11,v12,KAPPA,datathresh,probtfgene,sizeflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The PROM algorithm predicts the growth phenotype and the flux response
% after transcriptional perturbation, given a metabolic and regulatory
% network.
% INPUTS
% Model is the metabolic model obtained from COBRA toolbox through readcbmodel command
%
% Gene expression data - rows - genes,columns - conditions; no need to normalize or impute
%
% Expressionid - an array of identifiers for each row/gene should be included
%
% regulatory network - format - cell array of regulators and matching target genes
% example   Regulator = {'RegA'; 'RegB' ; 'RegC'};  Targets =
% {'GeneA';'GeneB';'GeneC'}
% note that the names or identifiers used in the regulatory data should
% match the names/ids given for gene expression data
%
%v11,v12 are obtained either from fastfva or fluxvariability command in
%cobra
% >> [v11,v12] = fastFVA(model);
%
%  sizeflag tells PROM if the regulatory network is large. it is 0 for large networks
%  and 1 for small networks ( less than 1000 interactions)
%
%OPTIONAL PARAMS
%
%- litevidence & prob_prior -> these should have the same length
%as the regulator/target arrays;
% high confidence interactions (not necessarily based on literature) should be flagged as one in litevidence
%array and the rest should be set as 0.
%
% Prob_prior should be set values between 0 and 1 for those interactions with litevidence. ( other values in
%the array would be ignored)
%
% KAPPA - determines strength of regulatory constraints - default value 1
% works for most systems
%
% subsets of tfs for which prom should be run ; default - for all tfs
%
% OUTPUT - the algorithm gives the growth rate (f) and flux response (v) after knock out of all
% the regulators in the regulatory model; status is the glpk solver status
%
% the status should be 5 for glpk; if its not then check solver error log
%
% lostxns gives the interactions that could not be quantified based on the
% threshold set for binarization
% the program would shoot a warning if the threshold chosen is bad. The
% default value (0.2 - 0.4) should work for most cases
%
% probtfgene gives the probabilities estimated for each interaction
%
%
% EXAMPLE
% load mtbpromdata
% f_ko = promv2(model,expression,expressionid,regulator,targets);
%[v11,v12] = fastFVA(model);
%[f,f_ko,v,v_ko,status1,lostxns,probtfgene] =  promv2(model,expression,expressionid,regulator,targets,[],[],[],v11,v12,[],[],[]);
% this is the tuberculosis model used in the PNAS paper
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT HANDLING
%===========================================================
if nargin == 5, litevidence = [];prob_prior = [];subsets = [];
elseif (nargin < 5) || (nargin == 6),
    error('Incorrect number of input arguments to %s',mfilename)
end
%===========================================================
%SOME BASIC INITIALIZATION
%===========================================================
disp('initializing data')

%% Parameter Initialization
% Tuning parameter KAPPA initialization
if (~exist('KAPPA','var')) || (isempty(KAPPA))
    KAPPA = 1;
end

% Presetting the DATATHRESHVAL parameter
if (~exist('DATATHRESHVAL','var')) || (isempty(DATATHRESHVAL))
    DATATHRESHVAL = 0.33; % KAPPA = 100;
end

%fprintf('params used - KAPPA: %d and DATATHRESH: %1.3f \n', KAPPA,DATATHRESHVAL)


%% TRN info initialization
regulated = targets;
[tfnames,b,m] = unique(regulator);


%% Metabolic Model variable initialization
weights = model.c;
S = model.S;
% ctype = repmat('S',size(model.b));
ctype = repmat('E',size(model.b));
dxdt = model.b;
param.msglev = 1;
lbff = model.lb;
ubff = model.ub;
% add a bit of wiggle room in flux bounds when upper = lower bound
lbff(lbff==ubff) = ubff(lbff == ubff) - 1E-6;

% number of reactions and metabolites in the model
numrxns = length(ubff);
nummets = size(S,1);

[rxnpos,genelist] = find(model.rxnGeneMat);
% u - reaction; v - genes
% finding rxn position

% if ~sizeflag
    if (~exist('v11','var')) || (isempty(v11))
        [v11,v12] = fluxVariability(model); 
        %[v11,v12] = fastFVA(model); 
    end
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find list of TFs to be considered in PROM knockout analysis
% i need to find the reactions that each gene has influence on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gene names from expression data
bnumsinexpsn = expressionid;

% genes with prior literature evidence
litevidence = logical(litevidence);

% Preparing to extract information on only a subset of TFs
if ~isempty(subsets)
    bnumstobekoed = subsets; % the TFs to be knocked out
    
    % extract the TF subnetwork containing the TF subset of interest
    regulated = targets(ismember(regulator,subsets));
    regulator = regulator(ismember(regulator,subsets));
       
else
    bnumstobekoed= tfnames;   % bnumstobekoed -  gives the geneids of the genes to be knocked out - by default it knocksout all the tfs in the model one by one
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scou = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Probabilities using a Global Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% remove_interactions1 = false(size(regulated));
disp('finding probabilities')
% cou3 = 1;

%% Quantile normalizing and imputing the expression data
data = expression;
data = knnimpute(data);
data = quantilenorm(data);  %its already normalized..
data1 = data;

%kappavals = [0,0.0001,0.001,0.05,0.1,0.25,0.33,0.5,0.75,1];
%datathreshvals = [0,0.01,0.05,0.1,0.2,0.25,0.33,0.4,0.5,0.75,1];
%datathresh = quantile(data(:),datathreshvals(scou));

%% Set the binarization threshold for the expression data. The variable
% 'datathreshval' gives the absolute expression vallue for the threshold;
% the variable DATATHRESHVAL gives the % quantile of the distribution

% if isempty(datathresh)
%     datathresh = quantile(data(:),DATATHRESHVAL);
% end
datathresh = quantile(data(:),DATATHRESHVAL);

% 'data' becomes the binarized version of the expression data here.
if datathresh < 0,
    data(data>=datathresh) = 1;
    data(data < datathresh) = 0;
else
    data(data < datathresh) = 0;
    data(data>=datathresh) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find probtfgene: Prob(Target = 1|TF = 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cou1 = 1;
cou = 1;
lost_xn = false(size(regulated));

if (~exist('probtfgene','var')) || (isempty(probtfgene))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try probtfgene = ones(size(regulated)); catch M1, probtfgene = 1; end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % iterates over each of the TF regulators in the TRN 
    for  i = 1:length(regulated)
        
        % finds which rows in the expression data correspond to the TF and
        % the target
        
        % k is the expression of the target gene being considered
        k = find(ismember(bnumsinexpsn,regulated(i)));
        % l is the TF expression
        l = find(ismember(bnumsinexpsn,regulator(i)));
        
        
        if ~isempty(k) & ~isempty(l)
            
            % te = target expression values across all samples
            targenExp = data1(k,:);
            
            % TF expression value across all samples
            tfExp = data1(l,:);
            
            % te = target binary expression values across all samples
            targenBin = data(k,:);
            
            % TF binary expression value across all samples
            tfBin = data(l,:);
            
            cou1 = cou1 + 1;
            
            % Kolmogorov-Smirnov test testing whether the two distributions
            % of target gene expressions (target expression where the TF is 'on'
            % and target expression when the TF is 'off') have significantly different 
            % expressions (does not necessarily assume same underlying
            % distribution/variance)
            try kstest2(targenExp(tfBin == 1),targenExp(tfBin== 0));
                
                if  (kstest2(targenExp(tfBin == 1),targenExp(tfBin== 0)) == 1),
                    
                    prob1 = sum(targenBin(tfBin == 0))/length(targenBin(tfBin==0));
                    probtfgene(i) = prob1;
                    cou = cou + 1;
                    %   this formula also gives the same answer  - (sum(~tec1(tec == 1))/length(tec1(tec==1))) * (sum(tec)/sum(~tec1))
                    
                else
                    
                    probtfgene(i) = 1;  % no effect
                    
                end
                
            catch ERRLG    % cant be estimated from microarray ; if it has strong evidence, i might consider setting this to zero later on
                probtfgene(i) = 1;
                lost_xn(i) = 1;
                
            end
            
        else
            probtfgene(i) = 1;
            lost_xn(i) = 1;
        end
    end
    
    probtfgene = probtfgene(:);
    
    toc
%    if (sum(lost_xn) > 0.75*length(probtfgene))   % check if there is problem with binarizartion.. usually there ll be some interactions that'd be missed, but if most of it (like >75% ) are...
%        % missed then its time to change the threshold
%        datathreshflag = 1;
%        disp('change binarization threshold')
%    else
        datathreshflag = 0;
%    end
    
    if ~isempty(litevidence)
        probtfgene(litevidence) = prob_prior(litevidence);  % u could set those interactions that u think have strong literature evidence to have predefined
    end                                                            % probabilities
    
    toc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  run PROM for each knockout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('running PROM')
count = 1;
thresh = 10^(-6);
mthresh = 10^(-3);
allgenes = [model.genes;tfnames];

% mapping the TRN target genes to the metabolic genes
[ir,posgenelist] = ismember(regulated,model.genes);

% standard FBA
lbf = lbff; ubf = ubff;
modeltmp.S = S;
modeltmp.c = weights;
modeltmp.b = dxdt;
modeltmp.lb = lbf;
modeltmp.ub = ubf;
modeltmp.csense = ctype;
modeltmp.stoic = S;
% [v,f_wt] = glpk(-weights,S,dxdt,lbf,ubf,ctype);
% f0 = -f_wt;
% solntmp = optimizeCbModelSM(modeltmp);
v = geometricFBASMprom(modeltmp,'epsilon',10^-8,'flexRel',10^-3);
% v = solntmp.x;
f0 = v(find(weights));

%%%%%%%
%% new additions
% Set up the optimization problem
    % max biomass function:
    % 1: S*v1 = 0
    % 3: v1 + delta+ >= -lb
    % 4: v1 - delta- <= ub
    % 5: -10^-6 <= delta+,delta- <= 0 
    %
    % 
% NOTE: lbg = lbff; ubg = ubff
lbg = model.lb; ubg = model.ub;

% add a bit of wiggle room in flux bounds when upper = lower bound
lbg(lbg==ubg) = ubg(lbg == ubg) - 1E-6;

% Creating an expanded stoichiometric matrix
a1 = [S,zeros(nummets,numrxns),zeros(nummets,numrxns)];
a2 = sparse([eye(numrxns), eye(numrxns),zeros(numrxns)]);
a3 = sparse([eye(numrxns), zeros(numrxns),-eye(numrxns)]);
A = [a1;a2;a3];

% creating an updated objective function
weights11 = [weights;zeros(2*numrxns,1)];
weights00 = [weights;zeros(2*numrxns,1)];

% creating updated upper- and lower- bound vectors
lb11 = [-1000*ones(numrxns,1);zeros(numrxns,1);zeros(numrxns,1)];
ub11 = [1000*ones(numrxns,1);zeros(numrxns,1);zeros(numrxns,1)];

% add a bit of wiggle room in flux bounds when upper = lower bound
lb11(lb11==ub11) = ub11(lb11 == ub11) - 1E-6;

% creating an updated model.b vector
dxdt0 = [zeros(nummets,1);lbg;ubg];

% ctype1 = [repmat('S',nummets,1);repmat('L',size(lbg,1),1);repmat('U',size(lbg,1),1)];
ctype1 = [repmat('E',nummets,1);repmat('G',size(lbg,1),1);repmat('L',size(lbg,1),1)];

% [x,f,origStat,extra] = glpk(c,A,b,lb,ub,csense,[],osense,params) is the
% standard optimizeCbModel formulation
% NOTE: v0 and f0 are the growth and flux profiles for the modified model
% created in this section
% [v0,f0] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);
% v = v0; %SM addition
%%%%%%


% flooring small flux values to zero
v(abs(v) < thresh) = 0;
v2 = nan(length(bnumstobekoed),numrxns);
% kapavals = [0,0.0001,0.001,0.05,0.1,0.25,0.33,0.5,1,5,10];
%kappa = kapavals(scou);
kappa = KAPPA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hw = waitbar(0);
% knock out the genes

vm = zeros(size(v));

% bnumstobekoed = [tfs_that_are_off;tf_that u want to knockout];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if ~sizeflag
    % disp('running fast fva for the entire network'); tic;
    %         [v11,v12] = fastFVA(model); toc;
    v11(abs(v11) < thresh) = 0;
    v12(abs(v12) < thresh) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% iterating over the TFs under consideration
for ci = 1:length(bnumstobekoed)
    
    %disp(ci)
    
    lbg = lbf; ubg = ubf;
    lb11 = [-1000*ones(numrxns,1);zeros(numrxns,1);zeros(numrxns,1)];
    ub11 = [1000*ones(numrxns,1);zeros(numrxns,1);zeros(numrxns,1)];

    % check if the TF being considered is also in the metabolic model
    if any(strcmpi(model.genes,bnumstobekoed(ci)))
        
        % temppos is the reactions that are mapped to the TF
        temppos = rxnpos(genelist == find(strcmp(model.genes,bnumstobekoed(ci))));
        
        % set upper and lower bounds of the reactions mapped to the TF to be 10^-6
        for jj = 1:length(temppos)
            if model.rev(temppos(jj))
                lbg(temppos) = -thresh;
                ubg(temppos) = thresh;
            else
                lbg(temppos) = -thresh;
            end
        end
        
    end
    
    
    % Redo FBA with the bounds adjusted for the rxns mapped to the TF
    % (similar to: [v,f_wt] = glpk(-weights,S,dxdt,lbf,ubf,ctype);)
    modeltmp.S = S;
    modeltmp.c = weights;
    modeltmp.b = dxdt;
    modeltmp.lb = lbg;
    modeltmp.ub = ubg;
    modeltmp.csense = ctype;
    modeltmp.stoic = S;
%     [v1,fk(ci)]  = glpk(-weights,S,dxdt,lbg,ubg,ctype);
%     solntmp = optimizeCbModelSM(modeltmp);
    v1 = geometricFBASMprom(modeltmp,'epsilon',10^-8,'flexRel',10^-3);
    fk(ci) = v1(find(weights));
%     fk(ci) = solntmp.f;
%     v1 = solntmp.x;
    
    % flooring small values to zero
    v1(abs(v1) < thresh) = 0; % need to check if this is actually used...
    
    % double check if the TF being considered is a TF in the TRN. 'tfnames'
    % is the set of unique TF regulators in the TRN.
    if any(ismember(tfnames,bnumstobekoed(ci)))      
        %  tfstate = false(size(tfnames));
%         tfstate(ismember(tfnames,bnumstobekoed(ci))) = 1;
%         
        
%         k = find(ismember(regulator,tfnames(tfstate)));
        k = find(ismember(regulator,bnumstobekoed(ci)));
        
        %    k(tempgeneprobs == 1) = '';
        
%       The target genes of the TF being considered, the probtf, and the
%       gene row and reaction row in the metabolic model
        tempgene = regulated(k);
        tempgeneprobs = probtfgene(k);
        tempgenepos = posgenelist(k);
        temprxnpos = rxnpos(ismember(genelist,tempgenepos));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% this section is for gene-protein-reaction relationship
        
        % x represents the metabolic gene state after a TF knockout,
        % assuming that all target genes are also knocked out with the TFko
        x = true(size(model.genes)); % corresponds to the rows of the genes in the metabolic model
        [isInModel,geneInd] = ismember(tempgene,model.genes);
%         x(geneInd) = false;
        x(geneInd(geneInd > 0)) = false; % x is false in the rows of the target genes (represents a target gene KO)
        
        constrainRxn = false(length(temprxnpos),1);
        % Figure out if any of the reaction states is changed
        % for the reactions mapped to the TF, find which are also knocked
        % out by model.rules
        for j = 1:length(temprxnpos)
            if (~eval(model.rules{temprxnpos(j)}))
                constrainRxn(j) = true; % records reactions that are perturbed
            end
        end
        % Constrain flux through the reactions associated with these genes
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Clear out the rows where ther are no metabolic genes mapped to the target
        % gene.
        tempgeneprobs(tempgenepos == 0)  = '';
        tempgene(tempgenepos == 0) = '';
        tempgenepos(tempgenepos == 0)  = '';
        
        
        
        % temprxnpos has the rxns that are going to  be affected by this tf
        % krxnpos are the rxns that will be affected by this target gene alone..
        % we loop around all the genes..
        
        % Iterating over the set of target genes that are in the metabolic model and mapped to the TF
        for l = 1:length(tempgenepos)
            
            % Check if the probtf is a reasonable number
            if ~isnan(tempgeneprobs(l))
                
                % Isolate the corresponding reactions mapped to the target
                % gene being considered
                krxnpos = ismember(temprxnpos,rxnpos(ismember(genelist,tempgenepos(l))));
                for m = 1:length(temprxnpos)
                    if krxnpos(m)
                        if constrainRxn(m)
                            % Consider only the reactions mapped to the
                            % target gene being considered that will be
                            % perturbed by the gene knockout, established
                            % by the GPR analysis section 
                            
                            % Consider only if probtf is between 0 and 1;
                            % if 0, result is like a direct knockout
                            if (tempgeneprobs(l) < 1)    % if its 1 no use in changing the bounds - might as well save time
                                if (tempgeneprobs(l) ~= 0)   % if its zero no point in estimating vm again - saves time.. but cant include in the above statement coz u have to change the bounds
                                    %if v(temprxnpos(m))
                                    
                                    if ~vm(temprxnpos(m))    % done to save time - if estimated already use it
                                        
%                                         % This section is skipped usually
%                                         % for large networks
%                                         if sizeflag
%                                             weights1 = weights; lbv = lbf; ubv = ubf;
%                                             grwthpos = find(weights == 1);
%                                             lbv(grwthpos) = v(grwthpos);
%                                             weights1(temprxnpos(m)) = -1;
%                                             [v11,fva1]  = glpk(-weights1,S,dxdt,lbv,ubv,ctype);
%                                             weights1(temprxnpos(m)) = 1;
%                                             [v12,fva2]  = glpk(-weights1,S,dxdt,lbv,ubv,ctype);
%                                             
%                                         end
                                        
                                        % The flux 'v' is from standard FBA
                                        % analysis. 'v11' is vmin from FVA,
                                        % 'v12' is vmax from FVA
                                        if  v(temprxnpos(m)) < 0
                                            % vm is the most extreme flux
                                            % value for that particular
                                            % reaction out of the FVA and
                                            % FBA flux values
                                            vm(temprxnpos(m)) = min([v11(temprxnpos(m)),v12(temprxnpos(m)),v(temprxnpos(m))]);
                                        elseif v(temprxnpos(m)) > 0
                                            vm(temprxnpos(m)) = max([v11(temprxnpos(m)),v12(temprxnpos(m)),v(temprxnpos(m))]);
                                        else
                                            vm(temprxnpos(m)) = max([abs(v11(temprxnpos(m))),abs(v12(temprxnpos(m))),abs(v(temprxnpos(m)))]);
                                        end
                                        
                                        
                                    end
                                    %end
                                end
                                
                                % Prepare the probtf adjustment of the flux bounds
                                xx =    vm(temprxnpos(m))*tempgeneprobs(l); % flux times probability
                                if  v(temprxnpos(m)) < 0
                                    
                                    % make sure that xx is less in
                                    % magnitude than the lbf and lbg
                                    % versions of the lower bounds
                                    tem = max([lbf(temprxnpos(m)),xx,lbg(temprxnpos(m))]);  %make sure we arent violating the original bounds; also get the lowest value if there were multiple modifications for the rxn
                                    
                                    % make sure that the absolute value of
                                    % the new flux bound is larger than
                                    % 10^-6
                                    lbg(temprxnpos(m)) = min([tem,-thresh]);   % prevents the solver from crashing
                                    
                                    ub11(1*numrxns + temprxnpos(m)) = 1000;
                                    weights11(1*numrxns + temprxnpos(m)) = ((-1)*kappa*abs(f0)/abs(vm(temprxnpos(m))));   % v0 f0 are the wild type values..
                                    vv = max([abs(vm(temprxnpos(m))),mthresh]);
                                    weights11(1*numrxns + temprxnpos(m)) = min([(-1)*kappa*abs(f0)/(abs(vv)), weights11(1*numrxns + temprxnpos(m)) ]);
                                elseif v(temprxnpos(m)) > 0
                                    tem = min([xx,ubf(temprxnpos(m)),ubg(temprxnpos(m))]);
                                    ubg(temprxnpos(m)) = max(tem,thresh);
                                    ub11(2*numrxns + temprxnpos(m)) = 1000;
                                    vv = max([abs(vm(temprxnpos(m))),mthresh]);
                                    weights11(2*numrxns + temprxnpos(m)) = min([(kappa*(-1)*abs(f0))/abs(vv), weights11(2*numrxns + temprxnpos(m)) ]);  % new weights based on kappa, normalized with growth rate
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    dxdt0 = [zeros(nummets,1);lbg;ubg];
    lb11(lb11==ub11) = ub11(lb11 == ub11) - 1E-6; % prevents solver from crashing

    modelmodtmp.S = A;
    modelmodtmp.c = weights11;
    modelmodtmp.b = dxdt0;
    modelmodtmp.lb = lb11;
    modelmodtmp.ub = ub11;
    modelmodtmp.csense = ctype1;
    modelmodtmp.stoic = S;
    %  FBA with the new bounds and weights
%     param.itlim = 1000000;
%     [v00(ci,:),f00(ci),status1(ci)] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);
    [vtmp stattmp]= geometricFBASMprom(modelmodtmp,'epsilon',10^-8,'flexRel',10^-3);
    if ~isempty(vtmp)
        v00(ci,:) = vtmp;
    else
        disp(['Fail in Line 577, ci = ' num2str(ci)])
        v00(ci,:) = nan(1,length(ub11));
    end
%     v00(ci,:) = vtmp;
    status1(ci) = stattmp;
%     v00(ci,:) = modsolntmp.x';
%     f00(ci) = modsolntmp.f;
%     status1(ci) = modsolntmp.origStat;
    
    f(ci) = v00(ci,find(weights));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coun  = 1;        lbh = lbg; ubh  = ubg;
    
    %% This section addresses cases where the updated FBA solution does not have a reaonable optimal solution
    while ((status1(ci) ~= 5) || (v00(ci,find(weights)) < 0))
        
        lbh(lbh ~= lbf) = lbh(lbh ~= lbf) - 1E-3;
        ubh(ubh ~= ubf) = ubh(ubh ~= ubf) + 1E-3;
        dxdt0 = [zeros(nummets,1);lbh;ubh];
        
        modelmodtmp.S = A;
        modelmodtmp.c = weights11;
        modelmodtmp.b = dxdt0;
        modelmodtmp.lb = lb11;
        modelmodtmp.ub = ub11;
        modelmodtmp.csense = ctype1;
        modelmodtmp.stoic = S;
%         [v00(ci,:),f00(ci),status1(ci)] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);
        [vtmp stattmp] = geometricFBASMprom(modelmodtmp,'epsilon',10^-8,'flexRel',10^-3);
        if ~isempty(vtmp)
            v00(ci,:) = vtmp';
        else
            disp(['Fail in Line 609, ci = ' num2str(ci)])
            v00(ci,:) = nan(1,length(ub11));
        end
        
        v00(ci,:) = vtmp;
        f00(ci) = v00(ci,find(weights));
        status1(ci) = stattmp;
        
        coun = coun + 1
        if (coun > 50),
            % if its impossible to estimate, then
            % check the unweighted g.r and the one with max weight prom -
            % if very less difference use it - any way warn the user about
            % the problem at the iteration number - ci;
            [v3,f3,status3] = glpk(-weights,S,dxdt,lbg,ubg,ctype);
            [v30,f30,status30] = glpk(-weights00,A,dxdt0,lb11,ub11,ctype1);
            %if abs((f3- v30(find(weights)))/abs(f3)) < 0.1
            f00(ci) = -f3;
            v00(ci,find(weights)) = -f3;
            %else
            disp(' problem in'); disp(ci);break;  % if that doesnt work,  display a warning
            %end
            
%             disp('check'); disp(ci); break;
            
        end
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f(ci) = v00(ci,find(weights));

    lbg_st(ci,:) = lbg;
    ubg_st(ci,:) = ubg;
    
    lb_st(ci,:) = lbff;
    ub_st(ci,:) = ubff;
    
    
    modeltmp2.S = S;
    modeltmp2.c = weights;
    modeltmp2.b = dxdt;
    modeltmp2.lb = lbg;
    modeltmp2.ub = ubg;
    modeltmp2.csense = ctype;
%     [v1,fk(ci)]  = glpk(-weights,S,dxdt,lbg,ubg,ctype);

%     [v2(ci,:),f1(ci),status] = glpk(-weights,S,dxdt,lbg,ubg,ctype);
    [vtmp status] = geometricFBASMprom(modelmodtmp2,'epsilon',10^-8,'flexRel',10^-3);
    
    if ~isempty(vtmp)
        v2(ci,:) = vtmp';
    end
    f1(ci) = vtmp(find(weights));
%     status = solntmp2.origStat;
    
%     f_ko(ci) = -f1(ci);
    f_ko(ci) = f1(ci);
    
%     ktime = toc;
%     waitpar = [num2str(ceil(ci/length(tfnames)*100)),'% complete. Time taken:',num2str(ceil(ktime)),' secs'];
    %waitbar(ci/length(tfnames),hw,waitpar);
    
    
    ff00(scou,ci) = v00(ci,find(weights));
    %disp(ff00(scou,ci))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%     this section is pretty much all ignored
%     if datathreshflag
%         if all(lost_xn(k))   % if none of the probabilities of a gene can be estimated, then ->
%             disp('Error: xns cant be estimated')
%             v00(ci,:) = NaN;
%             f1(ci) = NaN;
%             v2(ci,:) = NaN;
%             f00(ci) = NaN;
%             %break;
%         end
%     end
% %%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear tempgenepos tempgeneprobs temprxnpos k
    %disp(bnumstobekoed(ci))
end

% f_ko = -f1';
v_ko = v2;
%
% f00_ko(:,scou) = v00(:,find(weights));
v00_ko = v00;
%
% f = f00_ko;
v = v00_ko;


lostxns(:,scou) = lost_xn;

%disp([f,f_ko])

end

