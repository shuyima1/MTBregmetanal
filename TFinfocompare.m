for i = 1:size(uTF,1)
    TFname = uTF{i};
    
    eval(['siz = size(' TFname 'intInfo,1);'])
    TFstats(i,2) = siz - 1;
    if siz > 1
        eval(['RelgKO = cell2mat(' TFname 'intInfo(2:end,4))/soln.f;'])
        TFstats(i,3) = sum(RelgKO < 0.01);
        TFstats(i,4) = hygecdf(TFstats(i,3),size(iNJ661m.genes,1),grRat01,siz-1);
    end
end