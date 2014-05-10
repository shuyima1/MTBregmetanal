uTF2 = unique(chipnetOperon(:,1));
for i = 1:size(uTF2,1);

    TFname = uTF2{i};


    % Finds the TF-specific interactions from the TRN
    x = strmatch(TFname,chipnetOperon(:,1));
    eval([TFname 'ints = chipnetOperon(x,:);'])
    % isequal(' TFname 'ints(:,2),unique(' TFname 'ints(:,2),'stable'))

    % Finds the Metabolic Target genes
    eval(['[c ia ib] = intersect(' TFname 'ints(:,2),BesteSautonsJ.genes);'])
    eval([TFname 'Mints = ' TFname 'ints(ia,:);'])

    % Finds the ProbTFs associated with the TF of interest
    eval([TFname 'Pints = probtfgene(x);'])

    % Finds the ProbTFs associated with the target metabolic genes
    eval([TFname 'PMints = ' TFname 'Pints(ia);'])

    % Finds the list of affected Metabolic Target Genes and their corresponding
    % indices in the iNJ661m.genes (affectedix)
    eval([TFname 'affected = ' TFname 'PMints < 1;'])
    eval([TFname 'affectedix = [ib(' TFname 'affected) ' TFname 'PMints(' TFname 'affected)];'])


    % Finds the reactions mapped to the metabolic target genes
    eval(['flag = ~isempty(' TFname 'affectedix);'])
    if flag
        eval(['blah = BesteSautonsJ.rxnGeneMat(:,' TFname 'affectedix(:,1));'])
        eval([TFname 'RxnInfoJScnetOP = cell(size(blah,2),2);'])
        for j = 1:size(blah,2)
            tmp = find(blah(:,j) > 0);
            eval([TFname 'RxnInfoJScnetOP{j,1} = BesteSautonsJ.rxns(tmp);'])
            eval([TFname 'RxnInfoJScnetOP{j,2} = BesteSautonsJ.grRules(tmp);'])
        end
        eval([TFname 'RxnInfoJScnetOP(:,3) = BesteSautonsJ.genes(' TFname 'affectedix(:,1));'])
        eval([TFname 'RxnInfoJScnetOP(:,4) = num2cell(' TFname 'affectedix(:,2));'])
        eval(['tmpix = cellfun(@(x) strmatch(x,BesteSautonsJ.genes),' TFname 'RxnInfoJScnetOP(:,3));'])
        eval([TFname 'RxnInfoJScnetOP(:,5) = num2cell(grRatioJS(tmpix));'])

        % The long version of RxnInfo...has each reaction on a single line
        eval([TFname 'RxnLongJScnetOP = {};'])
        eval(['siz = size(' TFname 'RxnInfoJScnetOP,1);'])
        for j=1:siz
            eval(['cl = class(' TFname 'RxnInfoJScnetOP{j,1});'])
            if isequal(cl,'cell');
                eval([TFname 'RxnLongJScnetOP = [' TFname 'RxnLongJScnetOP;' TFname 'RxnInfoJScnetOP{j,1} repmat(' TFname 'RxnInfoJScnetOP(j,3),size(' TFname 'RxnInfoJScnetOP{j,1},1),1) repmat(' TFname 'RxnInfoJScnetOP(j,4),size(' TFname 'RxnInfoJScnetOP{j,1},1),1) repmat(' TFname 'RxnInfoJScnetOP(j,5),size(' TFname 'RxnInfoJScnetOP{j,1},1),1)];'])
            else
                eval([TFname 'RxnLongJScnetOP = [' TFname 'RxnLongJScnetOP;' TFname 'RxnInfoJScnetOP(j,[1 3 4 5])];'])
            end
        end
    end
    % Finds the Interaction Info: Regulator Target ProbTF grRateKO GriffinP Sassetti
    eval(['GrifmetIX = cellfun(@(x) strmatch(x,GriffinData(:,1),''exact''),' TFname 'Mints(:,2));'])
    eval([TFname 'intInfoJScnetOP = [' TFname 'Mints num2cell([' TFname 'PMints grRatioJS(ib)]) GriffinData(GrifmetIX,7:8)];'])
    eval([TFname 'intInfoJScnetOP = [{''Regulator'' ''Target'' ''ProbTF'' ''grRatioKO'' ''GriffinP'' ''Sassetti''};' TFname 'intInfoJScnetOP];'])

    eval(['save ' TFname 'AnalysisJScnetOP ' TFname '*'])
    
    fid = fopen([TFname 'intInfoJScnetOP2.txt'],'w');
    eval(['fprintf(fid, ''%s\t%s\t%s\t%s\t%s\t%s\n'',' TFname 'intInfoJScnetOP{1,:});'])
    
    eval(['nrows = size(' TFname 'intInfoJScnetOP,1);'])
    
    for j = 2:nrows
        eval(['fprintf(fid, ''%s\t%s\t%d\t%d\t%d\t%s\n'',' TFname 'intInfoJScnetOP{j,:});'])
    end
    fclose(fid);
end