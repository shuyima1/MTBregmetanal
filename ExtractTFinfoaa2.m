for i = 1:size(uTF,1);

    TFname = uTF{i};


    % Finds the TF-specific interactions from the TRN
    x = strmatch(TFname,BIGtf_05_2013_Rv2621(:,1));
    eval([TFname 'ints = BIGtf_05_2013_Rv2621(x,:);'])
    % isequal(' TFname 'ints(:,2),unique(' TFname 'ints(:,2),'stable'))

    % Finds the Metabolic Target genes
    eval(['[c ia ib] = intersect(' TFname 'ints(:,2),Beste7H9aa.genes);'])
    eval([TFname 'Mints = ' TFname 'ints(ia,:);'])

    % Finds the ProbTFs associated with the TF of interest
    eval([TFname 'Pints = probtfgeneB7aa(x);'])

    % Finds the ProbTFs associated with the target metabolic genes
    eval([TFname 'PMints = ' TFname 'Pints(ia);'])

    % Finds the list of affected Metabolic Target Genes and their corresponding
    % indices in the iNJ661m.genes (affectedix)
    eval([TFname 'affected = ' TFname 'PMints < 1;'])
    eval([TFname 'affectedix = [ib(' TFname 'affected) ' TFname 'PMints(' TFname 'affected)];'])


    % Finds the reactions mapped to the metabolic target genes
    eval(['flag = ~isempty(' TFname 'affectedix);'])
    if flag
        eval(['blah = Beste7H9aa.rxnGeneMat(:,' TFname 'affectedix(:,1));'])
        eval([TFname 'RxnInfoaa = cell(size(blah,2),2);'])
        for j = 1:size(blah,2)
            tmp = find(blah(:,j) > 0);
            eval([TFname 'RxnInfoaa{j,1} = Beste7H9aa.rxns(tmp);'])
            eval([TFname 'RxnInfoaa{j,2} = Beste7H9aa.grRules(tmp);'])
        end
        eval([TFname 'RxnInfoaa(:,3) = Beste7H9aa.genes(' TFname 'affectedix(:,1));'])
        eval([TFname 'RxnInfoaa(:,4) = num2cell(' TFname 'affectedix(:,2));'])
        eval(['tmpix = cellfun(@(x) strmatch(x,Beste7H9aa.genes),' TFname 'RxnInfoaa(:,3));'])
        eval([TFname 'RxnInfoaa(:,5) = num2cell(grRatioB7aa(tmpix));'])

        % The long version of RxnInfo...has each reaction on a single line
        eval([TFname 'RxnLongaa = {};'])
        eval(['siz = size(' TFname 'RxnInfoaa,1);'])
        for j=1:siz
            eval(['cl = class(' TFname 'RxnInfoaa{j,1});'])
            if isequal(cl,'cell');
                eval([TFname 'RxnLongaa = [' TFname 'RxnLongaa;' TFname 'RxnInfoaa{j,1} repmat(' TFname 'RxnInfoaa(j,3),size(' TFname 'RxnInfoaa{j,1},1),1) repmat(' TFname 'RxnInfoaa(j,4),size(' TFname 'RxnInfoaa{j,1},1),1) repmat(' TFname 'RxnInfoaa(j,5),size(' TFname 'RxnInfoaa{j,1},1),1)];'])
            else
                eval([TFname 'RxnLongaa = [' TFname 'RxnLongaa;' TFname 'RxnInfoaa(j,[1 3 4 5])];'])
            end
        end
    end
    % Finds the Interaction Info: Regulator Target ProbTF grRateKO GriffinP Sassetti
    eval(['GrifmetIX = cellfun(@(x) strmatch(x,GriffinData(:,1),''exact''),' TFname 'Mints(:,2));'])
    eval([TFname 'intInfoaa = [' TFname 'Mints num2cell([' TFname 'PMints grRatioB7aa(ib)]) GriffinData(GrifmetIX,7:8)];'])
    eval([TFname 'intInfoaa = [{''Regulator'' ''Target'' ''ProbTF'' ''grRatioKO'' ''GriffinP'' ''Sassetti''};' TFname 'intInfoaa];'])

    eval(['save ' TFname 'Analysisaa ' TFname '*'])
    
    fid = fopen([TFname 'intInfo2.txt'],'w');
    eval(['fprintf(fid, ''%s\t%s\t%s\t%s\t%s\t%s\n'',' TFname 'intInfoaa{1,:});'])
    
    eval(['nrows = size(' TFname 'intInfoaa,1);'])
    
    for j = 2:nrows
        eval(['fprintf(fid, ''%s\t%s\t%d\t%d\t%d\t%s\n'',' TFname 'intInfoaa{j,:});'])
    end
    fclose(fid);
end