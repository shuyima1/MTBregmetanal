uTF = unique(metregulator);

for i = 1:size(uTF,1);

    TFname = uTF{i};


    % Finds the TF-specific interactions from the TRN
%     x = strmatch(TFname,metregulator(:,1));
    eval([TFname 'ints = [metregulator(strcmp(TFname,metregulator(:,1))) metregulated(strcmp(TFname,metregulator(:,1)))];'])
    % isequal(' TFname 'ints(:,2),unique(' TFname 'ints(:,2),'stable'))

    % Finds the Metabolic Target genes
    eval(['[ia ia ib] = intersect(' TFname 'ints(:,2),Beste7H9Jaa2.genes);'])
    eval([TFname 'Mints = ' TFname 'ints(ia,:);'])

    % Finds the ProbTFs associated with the TF of interest
    eval([TFname 'Pints = probtfgenebootstrapstats(strcmp(TFname,metregulator(:,1)),:);'])

    % Finds the ProbTFs associated with the target metabolic genes
    eval([TFname 'PMints = ' TFname 'Pints(ia,:);'])

    % Finds the list of affected Metabolic Target Genes and their corresponding
    % indices in the iNJ661m.genes (affectedix)
    eval([TFname 'affected = ' TFname 'PMints(:,1) < 1;'])
    eval([TFname 'affectedix = [ib(' TFname 'affected) ' TFname 'PMints(' TFname 'affected,:)];'])


    % Finds the reactions mapped to the metabolic target genes
    eval(['flag = ~isempty(' TFname 'affectedix);'])
    if flag
        eval(['blah = Beste7H9Jaa2.rxnGeneMat(:,' TFname 'affectedix(:,1));'])
        eval([TFname 'RxnInfoBootstrap = cell(size(blah,2),2);'])
        for j = 1:size(blah,2)
            tmp = find(blah(:,j) > 0);
            eval([TFname 'RxnInfoBootstrap{j,1} = Beste7H9Jaa2.rxns(tmp);'])
            eval([TFname 'RxnInfoBootstrap{j,2} = Beste7H9Jaa2.grRules(tmp);'])
        end
        eval([TFname 'RxnInfoBootstrap(:,3) = Beste7H9Jaa2.genes(' TFname 'affectedix(:,1));'])
        eval([TFname 'RxnInfoBootstrap(:,4:6) = num2cell(' TFname 'affectedix(:,2:4));'])
        eval(['tmpix = cellfun(@(x) strmatch(x,Beste7H9Jaa2.genes),' TFname 'RxnInfoBootstrap(:,3));'])
        eval([TFname 'RxnInfoBootstrap(:,7) = num2cell(grRatioJaa2(tmpix));'])

        % The long version of RxnInfo...has each reaction on a single line
        eval([TFname 'RxnLongBootstrap = {};'])
        eval(['siz = size(' TFname 'RxnInfoBootstrap,1);'])
        for j=1:siz
            eval(['cl = class(' TFname 'RxnInfoBootstrap{j,1});'])
            if isequal(cl,'cell');
                eval([TFname 'RxnLongBootstrap = [' TFname 'RxnLongBootstrap;' TFname 'RxnInfoBootstrap{j,1} repmat(' TFname 'RxnInfoBootstrap(j,3),size(' TFname 'RxnInfoBootstrap{j,1},1),1) repmat(' TFname 'RxnInfoBootstrap(j,4:6),size(' TFname 'RxnInfoBootstrap{j,1},1),1) repmat(' TFname 'RxnInfoBootstrap(j,7),size(' TFname 'RxnInfoBootstrap{j,1},1),1)];'])
            else
                eval([TFname 'RxnLongBootstrap = [' TFname 'RxnLongBootstrap;' TFname 'RxnInfoBootstrap(j,[1 3 4 5 6 7])];'])
            end
        end
    end
    % Finds the Interaction Info: Regulator Target ProbTF grRateKO GriffinP Sassetti
    eval(['GrifmetIX = cellfun(@(x) strmatch(x,GriffinData(:,1),''exact''),' TFname 'Mints(:,2));'])
    eval([TFname 'intInfoBootstrap = [' TFname 'Mints num2cell([' TFname 'PMints grRatioJaa2(ib)]) GriffinData(GrifmetIX,7:8)];'])
    eval([TFname 'intInfoBootstrap = [{''Regulator'' ''Target'' ''ProbTF''  ''stdevProbTF''  ''CoefVarProbTF'' ''grRatioKO'' ''GriffinP'' ''Sassetti''};' TFname 'intInfoBootstrap];'])

    eval(['save ' TFname 'AnalysisBootstrap ' TFname '*'])
    
    fid = fopen([TFname 'intInfoBootstrap.txt'],'w');
    eval(['fprintf(fid, ''%s\t%s\t%s\t%s\t%s\t%s\n'',' TFname 'intInfoBootstrap{1,:});'])
    
    eval(['nrows = size(' TFname 'intInfoBootstrap,1);'])
    
    for j = 2:nrows
        eval(['fprintf(fid, ''%s\t%s\t%d\t%d\t%d\t%s\n'',' TFname 'intInfoBootstrap{j,:});'])
    end
    fclose(fid);
end