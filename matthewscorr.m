function [MCC] = matthewscorr(TP,TN,FP,FN)
% Calculates the Matthews Correlation Coefficient given input numbers (or
% vectors of identical size) of the TP, TN, FP, and FN

    if numel(TP) > 1
        MCC = ((TP.*TN)-(FP.*FN))./sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));
    elseif numel(TP) == 1
        MCC = ((TP*TN)-(FP*FN))./sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    else
        MCC = nan;
    end
end