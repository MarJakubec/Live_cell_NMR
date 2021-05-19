clear all 
close all 

load('Res_inten.mat')

% Convert time to hours and sort cell by it
time_tab = res_all(:,4);
time_f = sort(time_tab);

for i = 1:size(res_all,1)
t11=datevec(datenum(res_all(i,4)));
t22=datevec(datenum(time_f(1)));
time_dif = etime(t11,t22);
time_h = time_dif/3600;
res_all{i,5} = time_h;
end
out = sortrows(res_all,5);


min = cell2mat(out(:,5));

un_ppm_con = [];

for i = 1:size(out,1)
  cont_all(i,:) = out(i,1:5);
  un_ppm_con = [un_ppm_con; cont_all{i,2}];
end

un_ppm_con = unique(round(un_ppm_con,2));
overlap_con(size(un_ppm_con,1),size(min,1)) = 0;
for i = 1:size(un_ppm_con,1)
    for j = 1:size(min,1)
        ind =un_ppm_con(i) ==round(cont_all{j,2},2);
        if sum(ind) == 1
        con_fin(i,j) = cont_all{j,1}(ind);
        elseif sum(ind) > 1
            con_fin(i,j) = sum(cont_all{j,1}(ind));
            overlap_con(i,j) = 1;
            fprintf('Overlapping peaks i=%d j=%d \n',i,j)
        else 
        con_fin(i,j) = NaN;
        end
        clear ind
    end
end

name ='Res_plot';
mkdir(name)
addpath(name)
for i = 1:size(un_ppm_con,1)
figure('units','normalized','outerposition',[0 0 0.6 0.6],'visible','off');
hold on
scatter(min',con_fin(i,:),'LineWidth',2)
x = min';
y = con_fin(i,:);
xs = x(~isnan(y));
ys = y(~isnan(y));
plot(xs,ys,'r','LineWidth',2)
if any(overlap_con(i,:) == 1)
    ind = overlap_con(i,:) ==1 ;
x2 = x(ind);
y2 = y(ind);
scatter(x2,y2,'g','*','LineWidth',2)
f=get(gca,'Children');
legend(f(1),'Overlap','Location','southeast')
end
xlabel('Time [h]')
ylabel('Intensity')
% if max(y)>1E6
ylim([-Inf Inf])
% else
%     ylim([-Inf 1E6])
% end
xlim([0 55])
title(sprintf('%.2f ppm',un_ppm_con(i)))
set(gca,'Fontsize',14)
print(sprintf('%s/%s/Peak number %d',pwd,name,i),'-dpng')
close all
end
rmpath(name)


