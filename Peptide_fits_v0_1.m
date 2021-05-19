clear all
close all 

load('Guided.mat')
 
files = dir('*\pdata\*\ascii-spec.txt');

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


time = cell2mat(out(:,4));

un_ppm_con = [];
to_plot = [];
for i = 1:size(out,1)
to_plot = [to_plot out{i,1}];
end


nam = 'Pep_fit';
mkdir(nam)
addpath(nam)

%Peaks to fit - loop

to_fit = [47:48]; %position of peptide sin guided analysis

dim = [0.60 0.5 0.3 0.3]; %dimension and position of text 
m= 0;
for i = to_fit

x = time';
y = to_plot(i,:);
%WKWKWK - 0.85 uM
%  ([RL0]-[RLeq])*e^(kon*[ligand]+koff)*t)+[RLeq]
% correct (90374247-65041107)*2.718281828459046^-((kon*(0.85E-6)+koff)*x)+65041107
%equ  y = (MAX-MIN)*e^-(kon*0.4+koff)*x+MIN
%convert to concentration 
MAX = max(y);
y = y/MAX;
MAX = max(y);
MIN = mean(y(x>20));


[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( sprintf('(%d-%d)*2.718281828459046^-(kobs*x)+%d',MAX,MIN,MIN), 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.95];
opts.MaxIter = 1000000000;
opts.Robust = 'Bisquare';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
% gof
% fitresult
% 

yfit = fitresult([0:1:50]);
h = figure('units','normalized','outerposition',[0 0 0.7 0.8],'visible','off');
hold on
scatter(x,y,'filled')
plot([0:1:50],yfit,'r','LineWidth',2)
ylim([-Inf Inf])
set(gca,'Fontsize',18)
grid on
xlabel('time [h]')
ylabel('% of integral at zero')
title(sprintf('%s',peak_nam{i}))

kobs = fitresult.kobs; %uM/hod
t1_2 = 0.693/kobs;
t_eq = 5*t1_2;
% normalize 
FX = gradient(yfit);
slope = FX(1);
str = {sprintf('k_{obs}=%.3f',kobs),sprintf('t_{0.5}=%.3f',t1_2),sprintf('t_{eq}=%.3f',t_eq),...
    sprintf('slope_0=%.3f',slope),sprintf('r^2=%.3f',gof.rsquare)};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',16);
print(sprintf('%s/%s/%d_%s.png',pwd,nam,i,peak_nam{i}),'-dpng')
close all 
m= m+1;
res_slop(m) = slope;
res_kobs(m) = kobs;
res_teg(m) = t_eq;
res_thalf(m) = t1_2;
res_rsquare(m) = gof.rsquare; 
res_name{m} = peak_nam{i}; 
res_ate(m) = mean(y(x>t_eq));
end

save('Pep_fit.mat','res_slop','res_kobs','res_teg','res_thalf','res_rsquare',...
'res_name','res_ate')    
rmpath(nam)
