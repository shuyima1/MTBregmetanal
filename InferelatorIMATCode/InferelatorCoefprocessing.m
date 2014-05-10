MetCoefTF = cell(size(MetCoefMat,1),2);

for i = 1:size(MetCoefMat,1);
    x = strmatch(MetCoefMat{i,1},yeastY6trn(:,2),'exact');
    [c ia ib] = intersect(CoefRegulators,yeastY6trn(x,1),'stable');
    MetCoefTF{i,1} = CoefRegulators(ia);
    MetCoefTF{i,2} = cell2mat(MetCoefMat(i,ia+1))';
end