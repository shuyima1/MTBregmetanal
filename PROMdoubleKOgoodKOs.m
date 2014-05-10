% doubleJaa2goodKOs = [];
% for i = 1:size(fdoubleJaa2,1)
%     for j = 1:size(fdoubleJaa2,2)
%         if fdoubleJaa2(i,j) < mean(fdoubleJaa2(i,:)) - std(fdoubleJaa2(i,:))
%             if fdoubleJaa2(i,j) < mean(fdoubleJaa2(:,j)) - std(fdoubleJaa2(:,j))
%                 doubleJaa2goodKOs = [doubleJaa2goodKOs;i j fdoubleJaa2(i,j) mean(fdoubleJaa2(i,:)) mean(fdoubleJaa2(:,j))];
%             end
%         end
%     end
% end
% 
% doubleJaa2OPgoodKOs = [];
% for i = 1:size(fdoubleBIG50OPJaa2,1)
%     for j = 1:size(fdoubleBIG50OPJaa2,2)
%         if fdoubleBIG50OPJaa2(i,j) < mean(fdoubleBIG50OPJaa2(i,:)) - std(fdoubleBIG50OPJaa2(i,:))
%             if fdoubleBIG50OPJaa2(i,j) < mean(fdoubleBIG50OPJaa2(:,j)) - std(fdoubleBIG50OPJaa2(:,j))
%                 doubleJaa2OPgoodKOs = [doubleJaa2OPgoodKOs;i j fdoubleBIG50OPJaa2(i,j) mean(fdoubleBIG50OPJaa2(i,:)) mean(fdoubleBIG50OPJaa2(:,j))];
%             end
%         end
%     end
% end

% Synthetic sick
doubleJaa2goodKOs = [];
for i = 1:size(fdoubleJaa2,1)
    for j = 1:size(fdoubleJaa2,2)
        if fdoubleJaa2(i,j) < grKOJaa2(i);
            if fdoubleJaa2(i,j) < fJaa2(j);
                doubleJaa2goodKOs = [doubleJaa2goodKOs;i j fdoubleJaa2(i,j) grKOJaa2(i) fJaa2(j)]; % grKOJaa2 is the metabolic gene single deletion absolute growth rate, fJaa2 is the single TF deletion absolute growth rate
            end
        end
    end
end
doubleJaa2goodKOs(:,6) = doubleJaa2goodKOs(:,4) - doubleJaa2goodKOs(:,3);
doubleJaa2goodKOs(:,7) = doubleJaa2goodKOs(:,5) - doubleJaa2goodKOs(:,3);
doubleJaa2goodKOs(:,8) = mean(doubleJaa2goodKOs(:,6:7),2);
[ixb ixb] = sort(doubleJaa2goodKOs(:,8),'descend');
doubleJaa2goodKOsort = doubleJaa2goodKOs(ixb,:);

% Synthetic rescues
doubleJaa2goodKOs2 = [];
for i = 1:size(fdoubleJaa2,1)
    for j = 1:size(fdoubleJaa2,2)
        if fdoubleJaa2(i,j) > grKOJaa2(i);
            if fdoubleJaa2(i,j) > fJaa2(j);
                doubleJaa2goodKOs2 = [doubleJaa2goodKOs2;i j fdoubleJaa2(i,j) grKOJaa2(i) fJaa2(j)];
            end
        end
    end
end
doubleJaa2goodKOs2(:,6) = doubleJaa2goodKOs2(:,4) - doubleJaa2goodKOs2(:,3);
doubleJaa2goodKOs2(:,7) = doubleJaa2goodKOs2(:,5) - doubleJaa2goodKOs2(:,3);
doubleJaa2goodKOs2(:,8) = mean(doubleJaa2goodKOs2(:,6:7),2);
[ixb ixb] = sort(doubleJaa2goodKOs2(:,8),'ascend');
doubleJaa2goodKO2sort = doubleJaa2goodKOs2(ixb,:);

% doubleJaa2OPgoodKOs = [];
% for i = 1:size(fdoubleBIG50OPJaa2,1)
%     for j = 1:size(fdoubleBIG50OPJaa2,2)
%         if fdoubleBIG50OPJaa2(i,j) < mean(fdoubleBIG50OPJaa2(i,:)) - std(fdoubleBIG50OPJaa2(i,:))
%             if fdoubleBIG50OPJaa2(i,j) < mean(fdoubleBIG50OPJaa2(:,j)) - std(fdoubleBIG50OPJaa2(:,j))
%                 doubleJaa2OPgoodKOs = [doubleJaa2OPgoodKOs;i j fdoubleBIG50OPJaa2(i,j) mean(fdoubleBIG50OPJaa2(i,:)) mean(fdoubleBIG50OPJaa2(:,j))];
%             end
%         end
%     end
% end
