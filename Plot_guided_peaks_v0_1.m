clear all 
close all 

load('Guided.mat')

% Convert time to hours and sort cell by it
time_tab = res_qui(:,3);
time_f = sort(time_tab);

for i = 1:size(res_qui,1)
t11=datevec(datenum(res_qui(i,3)));
t22=datevec(datenum(time_f(1)));
time_dif = etime(t11,t22);
time_h = time_dif/3600;
res_qui{i,4} = time_h;
end
out = sortrows(res_qui,4);


min = cell2mat(out(:,4));

un_ppm_con = [];
to_plot = [];
for i = 1:size(out,1)
to_plot = [to_plot out{i,1}];
end
% %Delete outlying measurments 1 and thirs from end
% to_plot(:,1) = [];
% to_plot(:,35) = [];
% min(1) = [];
% min(35) = [];



% 
% 
% un_ppm_con = unique(round(un_ppm_con,2));
% overlap_con(size(un_ppm_con,1),size(min,1)) = 0;
% for i = 1:size(un_ppm_con,1)
%     for j = 1:size(min,1)
%         ind =un_ppm_con(i) ==round(cont_all{j,2},2);
%         if sum(ind) == 1
%         con_fin(i,j) = cont_all{j,1}(ind);
%         elseif sum(ind) > 1
%             con_fin(i,j) = sum(cont_all{j,1}(ind));
%             overlap_con(i,j) = 1;
%             fprintf('Overlapping peaks i=%d j=%d \n',i,j)
%         else 
%         con_fin(i,j) = NaN;
%         end
%         clear ind
%     end
% end

name ='Gui_plot';
mkdir(name)
addpath(name)
for i = 1:size(peak_nam,1)
figure('units','normalized','outerposition',[0 0 0.6 0.6],'visible','off');
hold on
x = min';
y = to_plot(i,:);
scatter(x,y,'LineWidth',2)
plot(x,y,'r','LineWidth',2)
xlabel('Time [h]')
ylabel('Integral')
if max(y)>1E6
ylim([-Inf Inf])
else
    ylim([-Inf 1E6])
end
xlim([0 55])
title(sprintf('%s %.2f-%.2f ppm',peak_nam{i},peak_int(i,1),peak_int(i,2)))
set(gca,'Fontsize',14)
print(sprintf('%s/%s/%d_%s.png',pwd,name,i,peak_nam{i}),'-dpng')
close all
end

% for i = 1:size(peak_nam,1)
% figure('units','normalized','outerposition',[0 0 0.6 0.6],'visible','off');
% hold on
% x = min(1:end-1)';
% y = to_plot(i,1:end-1);
% scatter(x,y,'LineWidth',2)
% plot(x,y,'r','LineWidth',2)
% xlabel('Time [h]')
% ylabel('Integral')
% ylim([-Inf Inf])
% xlim([0 45.2])
% title(sprintf('%s',peak_nam{i}))
% set(gca,'Fontsize',14)
% print(sprintf('%s/%s/%s_mod.png',pwd,name,peak_nam{i}),'-dpng')
% close all
% end
rmpath(name)
save('Guid_res.mat','to_plot','peak_nam','peak_int','min')

