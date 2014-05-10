PROMoldResults

thresh=[0:0.05:1]';
PROMoldRatio = cell2mat(PROMoldResults(2:end,1))/cell2mat(PROMoldResults(end,1));
PROMsassetti = PROMoldResults(2:end,2);
PROMgriffin = cell2mat(PROMoldResults(2:end,5));

PROMoldPerfsassetti = zeros(numel(thresh),3);
for i = 1:numel(thresh)
    PROMoldPerfsassetti(i,1) = sum(PROMoldRatio(strcmp('Essential',PROMsassetti)) < thresh(i))/sum(strcmp('Essential',PROMsassetti));
    PROMoldPerfsassetti(i,2) = sum(PROMoldRatio(strcmp('non-essential',PROMsassetti)) >= thresh(i))/sum(strcmp('non-essential',PROMsassetti));
end
PROMoldPerfsassetti(:,3) = PROMoldPerfsassetti(:,1)/2+PROMoldPerfsassetti(:,2)/2;


PROMoldSensgriffin = zeros(numel(thresh),3);
PROMoldSpecgriffin = zeros(numel(thresh),3);
truerange = [0.05 0.1 0.15]';
for i = 1:numel(thresh)
    for j = 1:numel(truerange)
        PROMoldSensgriffin(i,j) = sum(PROMoldRatio(PROMgriffin < 0.1) < thresh(i))/sum(PROMgriffin < 0.1);
        PROMoldSpecgriffin(i,j) = sum(PROMoldRatio(PROMgriffin >= 0.1) >= thresh(i))/sum(PROMgriffin >= 0.1);
    end
end

PROMoldPerfgriffin = PROMoldSensgriffin/2 + PROMoldSpecgriffin/2;