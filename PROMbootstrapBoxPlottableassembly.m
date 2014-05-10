load Griffin/GriffinData.mat
load('PROMbootstrapresults.mat')
load ../Beste7H9modJ BesteGriffinJ Beste7H9Jaa2
load mtbPROMinputs2 BIGt*

initCobraToolbox
sB7Jaa2 = optimizeCbModel(Beste7H9Jaa2);
sBG = optimizeCbModel(BesteGriffinJ);
fbootratioG = fbootG/sBG.f;
fbootstrapratio = fbootstrap/sB7Jaa2.f;


[c2 ia2 ib2] = intersect(unique(BIGtf_05_2013_Rv2621(:,1)),unique(metregulator));
[c ia ib] = intersect(unique(BIGtf_05_2013_Rv2621(:,1)),GriffinData(:,1));
c(:,2:3) = GriffinData(ib,7:8); %other columns of c are copied over from Excel.

fTFstats = {'TF','Griffin','Sassetti','NumRxns','Jaa2grRatio','GriffingrRatio'};
fid = fopen('PROMBIG50TFstats.txt','w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n',fTFstats{1,:});
for i = 1:50;fprintf(fid,'%s\t%d\t%s\t%d\t%d\t%d\n',c{i,:});end
fclose(fid);

fbootr2 = ones(500,50);
fbootr2(:,ia2) = fbootstrapratio';
fbootratioresults = {'TF','Griffin','Sassetti','BootstrapgrRatio'};
for i = 1:50;fbootratioresults = [fbootratioresults; repmat(c(i,:),500,1) num2cell(fbootr2(:,i))];end
fid = fopen('fbootstrapgrowthratiosJaa2.txt','w');
fprintf(fid,'%s\t%s\t%s\t%s\n',fbootratioresults{1,:});
for i = 2:size(fbootratioresults,1);fprintf(fid,'%s\t%d\t%s\t%d\n',fbootratioresults{i,:});end
fclose(fid);

fbootrG=ones(500,50);
fbootrG(:,ia2) = fbootratioG';
fbootratioresultsG = {'TF','BootstrapgrRatioGriffin'};
for i = 1:50;fbootratioresultsG = [fbootratioresultsG; repmat(c(i,1),500,1) num2cell(fbootrG(:,i))];end
fid = fopen('fbootstrapgrowthratiosGriffin.txt','w');
fprintf(fid,'%s\t%s\n',fbootratioresultsG{1,:});
for i = 2:size(fbootratioresultsG,1);fprintf(fid,'%s\t%d\n',fbootratioresultsG{i,:});end
fclose(fid);