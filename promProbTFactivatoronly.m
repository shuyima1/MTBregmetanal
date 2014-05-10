function [probtfgene, data, data1] =  promProbTFactivatoronly(model,expression,expressionid,regulator,targets,litevidence,prob_prior,datathresh,probtfgene)

%% INPUT HANDLING
%===========================================================
if nargin == 4, litevidence = [];prob_prior = [];
elseif (nargin < 5) || (nargin == 6),
    error('Incorrect number of input arguments to %s',mfilename)
end
%===========================================================
%SOME BASIC INITIALIZATION
%===========================================================
disp('initializing data')
regulated = targets;
% [tfnames] = unique(regulator);

% if (~exist('KAPPA','var')) || (isempty(KAPPA))
%     KAPPA = 1;
% end

if (~exist('DATATHRESHVAL','var')) || (isempty(DATATHRESHVAL))
    DATATHRESHVAL = 0.33; % KAPPA = 100;
end
%fprintf('params used - KAPPA: %d and DATATHRESH: %1.3f \n', KAPPA,DATATHRESHVAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

litevidence = logical(litevidence);
% if ~isempty(subsets)
%     bnumstobekoed = subsets;
%     regulated = targets(ismember(regulator,subsets));
%     regulator = regulator(ismember(regulator,subsets));
% else
%     bnumstobekoed= tfnames;   % bnumstobekoed -  gives the geneids of the genes to be knocked out - by default it knocksout all the tfs in the model one by one
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% new additions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Probabilities using a Global Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
lost_xn = false(size(regulated));
% remove_interactions1 = false(size(regulated));
disp('finding probabilities')
% cou = 1;cou1 = 1; %cou3 = 1;

data = expression;
data = knnimpute(data);
data = quantilenorm(data);  %its already normalized..

data1 = data;
%kappavals = [0,0.0001,0.001,0.05,0.1,0.25,0.33,0.5,0.75,1];
%datathreshvals = [0,0.01,0.05,0.1,0.2,0.25,0.33,0.4,0.5,0.75,1];
%datathresh = quantile(data(:),datathreshvals(scou));
if isempty(datathresh)
    datathresh = quantile(data(:),DATATHRESHVAL);
end


% datathresh = quantile(data(:),DATATHRESHVAL);

% disp(['datathresh=' num2str(datathresh)])

if datathresh < 0,
    data(data1>=datathresh) = 1;
    data(data1 < datathresh) = 0;
else
    data(data1 < datathresh) = 0;
    data(data1>=datathresh) = 1;
end

% data = databin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% numerator = ones(size(regulated));
% denominator = ones(size(regulated));
nboots = 500;

mettarints=cellfun(@(x) any(strcmp(x,model.genes)),regulated);

metregulated = regulated(mettarints);
metregulator = regulator(mettarints);

if (~exist('probtfgene','var')) || (isempty(probtfgene))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     try probtfgene = ones(size(regulated));
    try probtfgene = ones(size(metregulated,1),nboots);
    catch M1
        probtfgene = 1; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = cellfun(@(x) find(strcmp(x,expressionid)),metregulated);
    l = cellfun(@(x) find(strcmp(x,expressionid)),metregulator);
    for n = 1:nboots
        bootsamples = randsample(size(data,2),size(data,2),true);
        datatmp = data1(:,bootsamples);
        databintmp = data(:,bootsamples);
        for  i = 1:length(metregulated)
%             k = find(strcmp(metregulated(i),expressionid));
%             l = find(strcmp(metregulator(i),expressionid));
            if ~isempty(k(i)) & ~isempty(l(i))
                tarexp = datatmp(k(i),:);

    %             tfe1 = data1(l,:);

                tarbin = databintmp(k(i),:);
                tfbin = databintmp(l(i),:);

    %             cou1 = cou1 + 1;
                try ksh = kstest2(tarexp(tfbin == 1),tarexp(tfbin== 0));

                    if  ksh

                        prob1 = sum(tarbin(tfbin == 0))/length(tarbin(tfbin==0));
                        probOn = sum(tarbin(tfbin == 1))/length(tarbin(tfbin==1));
    %                     numerator(i) = sum(tarbin(tfbin == 0));
    %                     denominator(i) = length(tarbin(tfbin==0));
                        if probOn > prob1
                            probtfgene(i,n) = prob1;
                        end
    %                     cou = cou + 1;
                        %   this formula also gives the same answer  - (sum(~tfbin(tarbin == 1))/length(tfbin(tarbin==1))) * (sum(tarbin)/sum(~tfbin))

    %                 else
    %                     
    %                     probtfgene(i) = 1;  % no effect

                    end

                catch ERRLG    % cant be estimated from microarray ; if it has strong evidence, i might consider setting this to zero later on
    %                 probtfgene(i) = 1;
                    lost_xn(i) = 1;

                end
            else
    %             probtfgene(i) = 1;
                lost_xn(i) = 1;
            end
        end
    end
%     probtfgene = probtfgene(:);
    
    toc
%    if (sum(lost_xn) > 0.75*length(probtfgene))   % check if there is problem with binarizartion.. usually there ll be some interactions that'd be missed, but if most of it (like >75% ) are...
%        % missed then its time to change the threshold
%        datathreshflag = 1;
%        disp('change binarization threshold')
%    else
%         datathreshflag = 0;
%    end
    
    if ~isempty(litevidence)
        probtfgene(litevidence) = prob_prior(litevidence);  % u could set those interactions that u think have strong literature evidence to have predefined
    end                                                            % probabilities
    
    toc
end