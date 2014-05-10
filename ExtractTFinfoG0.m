for i = 1:size(uTF,1);

    TFname = uTF{i};


    % Finds the TF-specific interactions from the TRN
    x = strmatch(TFname,BIGtf_05_2013_Rv2621(:,1));
    eval([TFname 'ints = BIGtf_05_2013_Rv2621(x,:);'])
    % isequal(' TFname 'ints(:,2),unique(' TFname 'ints(:,2),'stable'))

    % Finds the Metabolic Target genes
    eval(['[c ia ib] = intersect(' TFname 'ints(:,2),BesteGriffin.genes);'])
    eval([TFname 'Mints = ' TFname 'ints(ia,:);'])

    % Finds the ProbTFs associated with the TF of interest
    eval([TFname 'Pints = probtfgeneG0(x);'])

    % Finds the ProbTFs associated with the target metabolic genes
    eval([TFname 'PMints = ' TFname 'Pints(ia);'])

    % Finds the list of affected Metabolic Target Genes and their corresponding
    % indices in the iNJ661m.genes (affectedix)
    eval([TFname 'affected = ' TFname 'PMints < 1;'])
    eval([TFname 'affectedix = [ib(' TFname 'affected) ' TFname 'PMints(' TFname 'affected)];'])


    % Finds the reactions mapped to the metabolic target genes
    eval(['flag = ~isempty(' TFname 'affectedix);'])
    if flag
        eval(['blah = BesteGriffin.rxnGeneMat(:,' TFname 'affectedix(:,1));'])
        eval([TFname 'RxnInfoG0 = cell(size(blah,2),2);'])
        for j = 1:size(blah,2)
            tmp = find(blah(:,j) > 0);
            eval([TFname 'RxnInfoG0{j,1} = BesteGriffin.rxns(tmp);'])
            eval([TFname 'RxnInfoG0{j,2} = BesteGriffin.grRules(tmp);'])
        end
        eval([TFname 'RxnInfoG0(:,3) = BesteGriffin.genes(' TFname 'affectedix(:,1));'])
        eval([TFname 'RxnInfoG0(:,4) = num2cell(' TFname 'affectedix(:,2));'])
        eval(['tmpix = cellfun(@(x) strmatch(x,BesteGriffin.genes),' TFname 'RxnInfoG0(:,3));'])
        eval([TFname 'RxnInfoG0(:,5) = num2cell(grRatioG(tmpix));'])

        % The long version of RxnInfo...has each reaction on a single line
        eval([TFname 'RxnLongG0 = {};'])
        eval(['siz = size(' TFname 'RxnInfoG0,1);'])
        for j=1:siz
            eval(['cl = class(' TFname 'RxnInfoG0{j,1});'])
            if isequal(cl,'cell');
                eval([TFname 'RxnLongG0 = [' TFname 'RxnLongG0;' TFname 'RxnInfoG0{j,1} repmat(' TFname 'RxnInfoG0(j,3),size(' TFname 'RxnInfoG0{j,1},1),1) repmat(' TFname 'RxnInfoG0(j,4),size(' TFname 'RxnInfoG0{j,1},1),1) repmat(' TFname 'RxnInfoG0(j,5),size(' TFname 'RxnInfoG0{j,1},1),1)];'])
            else
                eval([TFname 'RxnLongG0 = [' TFname 'RxnLongG0;' TFname 'RxnInfoG0(j,[1 3 4 5])];'])
            end
        end
    end
    % Finds the Interaction Info: Regulator Target ProbTF grRateKO GriffinP Sassetti
    eval(['GrifmetIX = cellfun(@(x) strmatch(x,GriffinData(:,1),''exact''),' TFname 'Mints(:,2));'])
    eval([TFname 'intInfoG0 = [' TFname 'Mints num2cell([' TFname 'PMints grRatioG(ib)]) GriffinData(GrifmetIX,7:8)];'])
    eval([TFname 'intInfoG0 = [{''Regulator'' ''Target'' ''ProbTF'' ''grRatioKO'' ''GriffinP'' ''Sassetti''};' TFname 'intInfoG0];'])

    eval(['save ' TFname 'AnalysisG0 ' TFname '*'])
    
    fid = fopen([TFname 'intInfo2.txt'],'w');
    eval(['fprintf(fid, ''%s\t%s\t%s\t%s\t%s\t%s\n'',' TFname 'intInfoG0{1,:});'])
    
    eval(['nrows = size(' TFname 'intInfoG0,1);'])
    
    for j = 2:nrows
        eval(['fprintf(fid, ''%s\t%s\t%d\t%d\t%d\t%s\n'',' TFname 'intInfoG0{j,:});'])
    end
    fclose(fid);
end