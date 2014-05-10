TPjaa2 = sum(BJaa2promstats.f(strcmp('essential',c(:,3)) | strcmp('growth-defect',c(:,3)))/sB7Jaa2.f < 0.85);
FNjaa2 = sum(BJaa2promstats.f(strcmp('essential',c(:,3)) | strcmp('growth-defect',c(:,3)))/sB7Jaa2.f >= 0.85);
TNjaa2 = sum(BJaa2promstats.f(strcmp('non-essential',c(:,3)))/sB7Jaa2.f >= 0.85);
FPjaa2 = sum(BJaa2promstats.f(strcmp('non-essential',c(:,3)))/sB7Jaa2.f < 0.85);
MCCjaa2 = ((TPjaa2*TNjaa2)-(FPjaa2*FNjaa2))/sqrt((TPjaa2+FPjaa2)*(TPjaa2+FNjaa2)*(TNjaa2+FPjaa2)*(TNjaa2+FNjaa2));

TPexpJaa2 = zeros(500,1);TNexpJaa2 = zeros(500,1); FPexpJaa2= zeros(500,1);FNexpJaa2 = zeros(500,1);
for i = 1:500;
TPexpJaa2(i,1)=sum(fnullexpperm(strcmp('essential',c(:,3)) | strcmp('growth-defect',c(:,3)),i)/sB7Jaa2.f < 0.85);
TNexpJaa2(i,1)=sum(fnullexpperm(strcmp('non-essential',c(:,3)),i)/sB7Jaa2.f >= 0.85);
FPexpJaa2(i,1)=sum(fnullexpperm(strcmp('non-essential',c(:,3)),i)/sB7Jaa2.f < 0.85);
FNexpJaa2(i,1)=sum(fnullexpperm(strcmp('essential',c(:,3)) | strcmp('growth-defect',c(:,3)),i)/sB7Jaa2.f >= 0.85);
end

TPnetJaa2 = zeros(500,1);TNnetJaa2 = zeros(500,1); FPnetJaa2= zeros(500,1);FNnetJaa2 = zeros(500,1);
for i = 1:500;
TPnetJaa2(i,1)=sum(fnullnetperm(strcmp('essential',c(:,3)) | strcmp('growth-defect',c(:,3)),i)/sB7Jaa2.f < 0.85);
TNnetJaa2(i,1)=sum(fnullnetperm(strcmp('non-essential',c(:,3)),i)/sB7Jaa2.f >= 0.85);
FPnetJaa2(i,1)=sum(fnullnetperm(strcmp('non-essential',c(:,3)),i)/sB7Jaa2.f < 0.85);
FNnetJaa2(i,1)=sum(fnullnetperm(strcmp('essential',c(:,3)) | strcmp('growth-defect',c(:,3)),i)/sB7Jaa2.f >= 0.85);
end

MCCexpjaa2 = ((TPexpJaa2.*TNexpJaa2)-(FPexpJaa2.*FNexpJaa2))./sqrt((TPexpJaa2+FPexpJaa2).*(TPexpJaa2+FNexpJaa2).*(TNexpJaa2+FPexpJaa2).*(TNexpJaa2+FNexpJaa2));

MCCnetjaa2 = ((TPnetJaa2.*TNnetJaa2)-(FPnetJaa2.*FNnetJaa2))./sqrt((TPnetJaa2+FPnetJaa2).*(TPnetJaa2+FNnetJaa2).*(TNnetJaa2+FPnetJaa2).*(TNnetJaa2+FNnetJaa2));


TPG = sum(BJGpromstats.f(cell2mat(c(:,2)) < 0.15)/sBG.f < 0.85);
FNG = sum(BJGpromstats.f(cell2mat(c(:,2)) < 0.15)/sBG.f >= 0.85);
TNG = sum(BJGpromstats.f(cell2mat(c(:,2)) >= 0.15)/sBG.f >= 0.85);
FPG = sum(BJGpromstats.f(cell2mat(c(:,2)) >= 0.15)/sBG.f < 0.85);

TPnetG = zeros(500,1);TNnetG = zeros(500,1); FPnetG= zeros(500,1);FNnetG = zeros(500,1);
for i = 1:500;
TPnetG(i,1)=sum(fnullnetpermG(cell2mat(c(:,2)) < 0.15,i)/sBG.f < 0.85);
TNnetG(i,1)=sum(fnullnetpermG(cell2mat(c(:,2)) >= 0.15,i)/sBG.f >= 0.85);
FPnetG(i,1)=sum(fnullnetpermG(cell2mat(c(:,2)) >= 0.15,i)/sBG.f < 0.85);
FNnetG(i,1)=sum(fnullnetpermG(cell2mat(c(:,2)) < 0.15,i)/sBG.f >= 0.85);
end

TPexpG = zeros(500,1);TNexpG = zeros(500,1); FPexpG= zeros(500,1);FNexpG = zeros(500,1);
for i = 1:500;
TPexpG(i,1)=sum(fnullexppermG(cell2mat(c(:,2)) < 0.15,i)/sBG.f < 0.85);
TNexpG(i,1)=sum(fnullexppermG(cell2mat(c(:,2)) >= 0.15,i)/sBG.f >= 0.85);
FPexpG(i,1)=sum(fnullexppermG(cell2mat(c(:,2)) >= 0.15,i)/sBG.f < 0.85);
FNexpG(i,1)=sum(fnullexppermG(cell2mat(c(:,2)) < 0.15,i)/sBG.f >= 0.85);
end