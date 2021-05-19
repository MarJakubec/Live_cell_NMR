clear all
close all 

% load('Res_inten.mat')
 % guided peaks
 % sucrose - 5.300 -5.336; 3.354-3.408
 % sukcinate 2.303 - 2.313; 
 % metabolites? 1.086-1.104; 1.148-1.158
 
 files = dir('*\pdata\*\ascii-spec.txt');
% files_names = dir('*\pdata\*\title');
% fid_files = dir('*\fid');


for i = 1:size(files,1)
addpath(files(i).folder)
out{i} = dlmread(files(i).name,',', 1, 0);
% Import identifier
name{i} = fileread('title');
rmpath(files(i).folder) 
end

peak_int = [2.3009 2.3126;
    1.8183 1.8240;
    2.0824 2.0963;
    2.5174 2.6045;
    3.3472 3.3779;
    7.757 7.8650;
    8.0887 8.1044;
    8.1044 8.1212;
    8.35447 8.36462;
    ...%lyzate_SPECIFIC
    0.9538 0.9650;
    0.9650 0.9767;
    1.3690 1.3800;
    5.6957 5.7203;
    7.4226 7.4512;
    8.0814 8.1012;
    8.1012 8.1202;
    8.2385 8.2476;
    0.8488 0.8585;
    0.8585 0.8699;
    0.9468 0.9540;
    0.9548 0.9650;
    1.7667 1.7790;
    1.9341 1.9417;
    ...%AMP16_lys_spec
    0.8479 0.8688;
    0.8688 0.8788;
    0.8788 0.9015;
    0.9015 0.9245;
    0.9296 0.9567;
    2.2367 2.2761;
    2.8950 2.9068;
    2.9068 2.9188;
    2.9188 2.9321;
    3.0790 3.0910;
    3.1359 3.1708;
    7.1010 7.1100;
    7.17328 7.1855;
    7.1855 7.19535;
    7.19535 7.2072;
    7.2164 7.2322;
    7.4363 7.4621;
    7.6247 7.6543;
    ...%AMP26_lys_specific
    1.5848 1.6471;
    2.9059 2.9466;
    3.1625 3.1710;
    ...%AMPicilin SPECIFIC? 
    2.0859 2.0951;
    3.2327 3.2451;
    ...%PEPTIDE SPECIFIC
    2.6152 2.6409;
    3.2581 3.2698];

peak_nam = {'Sukcinate'; 
    'Metabolite singlet';
    'Metabolite singlet_AMP21';
    'Metabolite triplet';
    'Metabolite triplet';
    'Metabolite singlet';
    'Metabolite singlet ';
    'Metabolite singlet';
    'Metabolite singlet';
    ...%LYZATE Specific 
    'Lys S';
    'Lys S';
    'Disruption of bilayer S';
    'Disruption of bilayer D';
    'Disruption of bilayer  D';
    'Disruption of bilayer  S';
    'Disruption of bilayer  S';
    'Disruption of bilayer  S';
    'Lys S';
    'Lys S';
    'Lys S';
    'Lys S';
    'Lys S';
    'Lys S';
    ...%AMP 16 lyz specific
    'AMP 16 lys D';
    'AMP 16 lys S';
    'AMP 16 lys D';
    'AMP 16 lys D';
    'AMP 16 lys D';
    'AMP 16 lys M';
    'AMP 16 lys S';
    'AMP 16 lys S';
    'AMP 16 lys S';
    'AMP 16 lys S'; 
    'AMP 16 lys T';
    'AMP 16 lys T'; 
    'AMP 16 lys M';
    'AMP 16 lys M';
    'AMP 16 lys M';
    'AMP 16 lys S';
    'AMP 16 lys D';
    'AMP 16 lys DD';
    ...%AMP 26 lys specific
    'AMP 26 lys P';
    'AMP 26_AMP 21 lys T';
    'AMP 26 lys S';
    ...%AMP specific
    'Ampicilin Met S';
    'Ampicilin Met S';
   ...%Peptide specific
   'AMP 95 S';
   'AMP 95 M'};

for i = 1:size(out,2)
mtest(:,1) = out{i}(:,4);
mtest(:,2) = out{i}(:,2);
% mtest2 = mtest(out{i}(:,4) > 0 & out{i}(:,4)<3.5,:);
% mtest = mtest(out{i}(:,4) > 5 & out{i}(:,4)<11,:);
% mtest = flipud(mtest);
% mtest2 = flipud(mtest2);
x = mtest(:,1); 
y = mtest(:,2);
% noise estimate
nos = mtest(mtest(:,1)>10.5,2);
nos_v = abs(min(nos))+ abs(max(nos));
nos_v = 1.5*nos_v;
% figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
for j = 1:size(peak_int,1)
%     siz = peak_int(j,2) - peak_int(j,1); 
    ind = x >= peak_int(j,1) & x <= peak_int(j,2);
    int(j,1) = trapz(y(ind,:));
%     ind2 =  x > peak_int(j,1)-siz & x<peak_int(j,2)+siz;
%     subplot(size(peak_int,1),1,j)
%     plot(x(ind2),y(ind2))
%     hold on
%     plot(x(ind),zeros(size(x(ind),1)),'r','LineWidth',1.5)
%     set(gca, 'XDir','reverse')
%     title(peak_nam{j})
%     set(gca,'Fontsize',12)
%     hold off  
end
to_scan = files(i).folder(end-12:end);
ind = regexp(to_scan,'\\');
th = sscanf(to_scan(ind(1):ind(2)),'%*[^0123456789]%d');
% sgtitle(name{i})
% print(sprintf('Guided_analysis %d.png',th),'-dpng')
close all


res_qui{i,1} = int;
res_qui{i,2} = name(i);
fid_fil = dir(sprintf('%s/fid',files(i).folder(1:end-7)));
res_qui{i,3} = fid_fil.date;
clear mtest int
end

save('Guided.mat','res_qui','peak_nam','peak_int')

clear all 
close all 
