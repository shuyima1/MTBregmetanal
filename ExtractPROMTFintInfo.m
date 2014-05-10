%% ExtractPROMTFintInfo.m
% Shuyi Ma 2014

% This code extracts the information for each TF associated with its target
% metabolic genes, the probability of the interaction, and the single gene
% deletion values of the target metabolic genes. Useful for interpreting
% which target metabolic genes are primarily responsible for any growth
% deficit predicted for the TF.

% Inputs:
% model: Metabolic Model
% regulators: (same as 'regulator' input variable from promv2). Cell array of length equal to the number of interactions in the TRN listing the regulator of each interaction in the TRN
% targets: (same as 'regulated' input from promv2). Cell array of length equal to the number of interactions in the TRN listing the target of each interaction in the TRN
% probtfgene: (probtfgene output from promv2).
% singleGeneDelGrowthRatio: single metabolic gene deletion growth ratio from COBRA 
% Experimental Essentiality: optional cell array with one column of genes
% and a second column of essentiality data for those genes. If unavailable,
% use [] as input.
% saveflag: 1 if save the output, 0 if no saving
% textoutflag: 1 if write interaction information to text file, 0 otherwise

% Outputs:
% TFintInfo: cell array containing the information of each of the
% interactions mapped to each TF (Regulator, Target, ProbTF, MetabolicGeneKO
% Ratio, Experimental Essentiality Metric
% TFrxnLong: Cell array containing reaction mapping information for each TF

function [TFintInfo, TFrxnLong] = ExtractPROMTFintInfo(model,regulators,targets,probtfgene,singleGeneDelGrowthRatio,ExperimentalEssentialityMetric,saveflag,textoutflag)


uTF = unique(regulators(:,1));

TFintInfo = cell(size(uTF,1),1);
TFrxnLong = cell(size(uTF,1),1);
% singleGeneDelGrowthRatio = singleGeneDeletion(model);
% [probtfgene probtfgene probtfgene probtfgene probtfgene probtfgene probtfgene] = promv2(model,RvFCexp,Rvgens(:,1),regulators,AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);

for i = 1:size(uTF,1);

    TFname = uTF{i};

    % Finds the TF-specific interactions from the TRN
    x = strcmp(TFname,regulators);
    
    TFints = [regulators(x) targets(x)];
    % isequal(' TFname 'ints(:,2),unique(' TFname 'ints(:,2),'stable'))

    % Finds the ProbTFs associated with the TF of interest
    intsProbTF = probtfgene(x);
    
    % Finds the Metabolic Target genes
    [ia, ia, ib] = intersect(TFints(:,2),model.genes);
    TFMints = TFints(ia,:);

    % Finds the ProbTFs associated with the target metabolic genes
    MetintsProbTF = intsProbTF(ia);

    % Finds the list of affected Metabolic Target Genes and their corresponding
    % indices in the metabolic model (TFaffectedix)
    TFaffectedprobinfo = [ib(MetintsProbTF < 1) MetintsProbTF(MetintsProbTF < 1)];

    % Finds the reactions mapped to the metabolic target genes
    if ~isempty(TFaffectedprobinfo)
        blah = model.rxnGeneMat(:,TFaffectedprobinfo(:,1));
        TFRxnInfo = cell(size(blah,2),2);
        for j = 1:size(blah,2)
            tmp = find(blah(:,j) > 0);
            TFRxnInfo{j,1} = model.rxns(tmp);
            TFRxnInfo{j,2} = model.grRules(tmp);
        end
        TFRxnInfo(:,3) = model.genes(TFaffectedprobinfo(:,1));
        TFRxnInfo(:,4) = num2cell(TFaffectedprobinfo(:,2));
        tmpix = cellfun(@(x) strmatch(x,model.genes),TFRxnInfo(:,3));
        TFRxnInfo(:,5) = num2cell(singleGeneDelGrowthRatio(tmpix));

        % The long version of RxnInfo...has each reaction on a single line
        for j=1:size(TFRxnInfo,1)
            if isequal(class(TFRxnInfo{j,1}),'cell');
                TFrxnLong{i} = [TFrxnLong{i};TFRxnInfo{j,1} repmat(TFRxnInfo(j,3),size(TFRxnInfo{j,1},1),1) repmat(TFRxnInfo(j,4),size(TFRxnInfo{j,1},1),1) repmat(TFRxnInfo(j,5),size(TFRxnInfo{j,1},1),1)];
            else
                TFrxnLong{i} = [TFrxnLong{i};TFRxnInfo(j,[1 3 4 5])];
            end
        end
    end
    
    if ~isempty(ExperimentalEssentialityMetric)
    % Finds the Interaction Info: Regulator Target ProbTF grRateKO GriffinP Sassetti
        EssentialitymetIX = cellfun(@(x) find(strcmp(x,ExperimentalEssentialityMetric(:,1))),TFMints(:,2));
    
        TFintInfo{i} = [{'Regulator' 'Target' 'ProbTF' 'grRatioKO' 'ExperimentalEssentiality'};TFMints num2cell([MetintsProbTF singleGeneDelGrowthRatio(ib)]) ExperimentalEssentialityMetric(EssentialitymetIX,2)];
    else
        TFintInfo{i} = [{'Regulator' 'Target' 'ProbTF' 'grRatioKO'};TFMints num2cell([MetintsProbTF singleGeneDelGrowthRatio(ib)])];
    end
    
    if saveflag
        eval([TFname 'intInfo = TFintInfo{i};'])
        eval([TFname 'RxnLong = TFrxnLong{i};'])
        eval(['save ' TFname 'AnalysisAltPeakNetOPJS ' TFname 'intInfo ' TFname 'RxnLong'])
    end
    
    if textoutflag
        fid = fopen([TFname 'intInfo.txt'],'w');
        fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n',TFintInfo{i}{1,:});
    
        for j = 2:size(TFintInfo{i},1)
            eval(['fprintf(fid, ''%s\t%s\t%d\t%d\t%d\t%s\n'',' TFname 'intInfo{i}{j,:});'])
        end
        fclose(fid);
    end
end

end