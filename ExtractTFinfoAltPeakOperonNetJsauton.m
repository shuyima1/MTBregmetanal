%load AltPeakNetwork AltPeakOperonNet
%load Beste7H9modJ BesteSautonsJ
%load PROMchipnetInputs RvFCexp Rvgens

initCobraToolbox
changeCobraSolver('glpk')
uTF = unique(AltPeakOperonNet(:,1));

grRatioJS = singleGeneDeletion(BesteSautonsJ);
[probtfgene probtfgene probtfgene probtfgene probtfgene probtfgene probtfgene] = promv2(BesteSautonsJ,RvFCexp,Rvgens(:,1),AltPeakOperonNet(:,1),AltPeakOperonNet(:,2),[],[],[],[],[],[],0,[],0);

for i = 1:size(uTF,1);

    TFname = uTF{i};


    % Finds the TF-specific interactions from the TRN
    x = strmatch(TFname,AltPeakOperonNet(:,1));
    eval([TFname 'ints = AltPeakOperonNet(x,:);'])
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
        eval([TFname 'RxnInfoAPJS = cell(size(blah,2),2);'])
        for j = 1:size(blah,2)
            tmp = find(blah(:,j) > 0);
            eval([TFname 'RxnInfoAPJS{j,1} = BesteSautonsJ.rxns(tmp);'])
            eval([TFname 'RxnInfoAPJS{j,2} = BesteSautonsJ.grRules(tmp);'])
        end
        eval([TFname 'RxnInfoAPJS(:,3) = BesteSautonsJ.genes(' TFname 'affectedix(:,1));'])
        eval([TFname 'RxnInfoAPJS(:,4) = num2cell(' TFname 'affectedix(:,2));'])
        eval(['tmpix = cellfun(@(x) strmatch(x,BesteSautonsJ.genes),' TFname 'RxnInfoAPJS(:,3));'])
        eval([TFname 'RxnInfoAPJS(:,5) = num2cell(grRatioJS(tmpix));'])

        % The long version of RxnInfo...has each reaction on a single line
        eval([TFname 'RxnLongAPJS = {};'])
        eval(['siz = size(' TFname 'RxnInfoAPJS,1);'])
        for j=1:siz
            eval(['cl = class(' TFname 'RxnInfoAPJS{j,1});'])
            if isequal(cl,'cell');
                eval([TFname 'RxnLongAPJS = [' TFname 'RxnLongAPJS;' TFname 'RxnInfoAPJS{j,1} repmat(' TFname 'RxnInfoAPJS(j,3),size(' TFname 'RxnInfoAPJS{j,1},1),1) repmat(' TFname 'RxnInfoAPJS(j,4),size(' TFname 'RxnInfoAPJS{j,1},1),1) repmat(' TFname 'RxnInfoAPJS(j,5),size(' TFname 'RxnInfoAPJS{j,1},1),1)];'])
            else
                eval([TFname 'RxnLongAPJS = [' TFname 'RxnLongAPJS;' TFname 'RxnInfoAPJS(j,[1 3 4 5])];'])
            end
        end
    end
    % Finds the Interaction Info: Regulator Target ProbTF grRateKO GriffinP Sassetti
    eval(['GrifmetIX = cellfun(@(x) strmatch(x,GriffinData(:,1),''exact''),' TFname 'Mints(:,2));'])
    eval([TFname 'intInfoAPJS = [' TFname 'Mints num2cell([' TFname 'PMints grRatioJS(ib)]) GriffinData(GrifmetIX,7:8)];'])
    eval([TFname 'intInfoAPJS = [{''Regulator'' ''Target'' ''ProbTF'' ''grRatioKO'' ''GriffinP'' ''Sassetti''};' TFname 'intInfoAPJS];'])

    eval(['save ' TFname 'AnalysisAltPeakNetOPJS ' TFname '*'])
    
    fid = fopen([TFname 'intInfoAltPeakNetOPJS.txt'],'w');
    eval(['fprintf(fid, ''%s\t%s\t%s\t%s\t%s\t%s\n'',' TFname 'intInfoAPJS{1,:});'])
    
    eval(['nrows = size(' TFname 'intInfoAPJS,1);'])
    
    for j = 2:nrows
        eval(['fprintf(fid, ''%s\t%s\t%d\t%d\t%d\t%s\n'',' TFname 'intInfoAPJS{j,:});'])
    end
    fclose(fid);
end