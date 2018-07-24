function my_plot()
close all;
clear; clc;
IniToolbox;

% ordinary_3d_results();
% planar_3d_results();
odinary_3d_results_sigma();
planar_3d_results_sigma();
% odinary_3d_results_time();
end
%%
function ordinary_3d_results()
load ordinary3Dresults
figure('outerposition',[0    260    1920    400])
yrange= [0.1 1.1];

subplot(1,2,1)
grid on
xdrawgraph(npts,yrange,method_list,'deleted_mean_r','Mean Rotation Error',...
    'Number of Points','Rotation Error (degrees)');
legend('hide')
ca_pos = get(gca,'OuterPosition');
ca_pos(4)=ca_pos(4)-0.1;
set(gca,'OuterPosition', ca_pos);

subplot(1,2,2)
grid on
xdrawgraph(npts,yrange,method_list,'deleted_mean_t','Mean Translation Error',...
    'Number of Points','Translation Error (%)')
%legend('hide')
ca_pos = get(gca,'OuterPosition');
ca_pos(4)=ca_pos(4)-0.1;
set(gca,'OuterPosition', ca_pos);


hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend,'Location','northoutside','Orientation','horizontal','position',[0.5 0.95 0 0]);

saveas(gcf, 'sim_experiment_ordinary', 'epsc');
end
%%
function planar_3d_results()
load planar3Dresults
figure('outerposition',[0    260    1920    400])
yrange= [0.1 2];

subplot(1,2,1)
grid on
xdrawgraph(npts,yrange,method_list,'deleted_mean_r','Mean Rotation Error',...
    'Number of Points','Rotation Error (degrees)');
legend('hide')
ca_pos = get(gca,'OuterPosition');
ca_pos(4)=ca_pos(4)-0.1;
set(gca,'OuterPosition', ca_pos);

subplot(1,2,2)
yrange= [0.1 1.2];
grid on
xdrawgraph(npts,yrange,method_list,'deleted_mean_t','Mean Translation Error',...
    'Number of Points','Translation Error (%)')
ca_pos = get(gca,'OuterPosition');
ca_pos(4)=ca_pos(4)-0.1;
set(gca,'OuterPosition', ca_pos);


hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend,'Location','northoutside','Orientation','horizontal','position',[0.5 0.95 0 0]);
legend('hide')

saveas(gcf, 'sim_experiment_planar' ,'epsc')
end
%%
function odinary_3d_results_sigma()
load ordinary3DresultsSigma
figure('outerposition',[0    260    1920    400])
yrange= [0 0.8];

subplot(1,2,1)
grid on
xdrawgraph(nls,yrange,method_list,'deleted_mean_r','Mean Rotation Error',...
    'Gaussian Image Noise (pixels)','Rotation Error (degrees)');
legend('hide')
ca_pos = get(gca,'OuterPosition');
ca_pos(4)=ca_pos(4)-0.1;
set(gca,'OuterPosition', ca_pos);

subplot(1,2,2)
grid on
xdrawgraph(nls,yrange,method_list,'deleted_mean_t','Mean Translation Error',...
    'Gaussian Image Noise (pixels)','Translation Error (%)');
ca_pos = get(gca,'OuterPosition');
ca_pos(4)=ca_pos(4)-0.1;
set(gca,'OuterPosition', ca_pos);


hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend,'Location','northoutside','Orientation','horizontal','position',[0.5 0.95 0 0]);

saveas(gcf, 'sim_experiment_ordinary_sigma' ,'epsc')
end
%%
function planar_3d_results_sigma()
load planar3DresultsSigma
figure('outerposition',[0    260    1920    400])
yrange= [0 1.3];
subplot(1,2,1)
grid on
xdrawgraph(nls,yrange,method_list,'deleted_mean_r','Mean Rotation Error',...
    'Gaussian Image Noise (pixels)','Rotation Error (degrees)');
legend('hide')
ca_pos = get(gca,'OuterPosition');
ca_pos(4)=ca_pos(4)-0.1;
set(gca,'OuterPosition', ca_pos);

yrange= [0 0.6];
subplot(1,2,2)
grid on
xdrawgraph(nls,yrange,method_list,'deleted_mean_t','Mean Translation Error',...
    'Gaussian Image Noise (pixels)','Translation Error (%)');
ca_pos = get(gca,'OuterPosition');
ca_pos(4)=ca_pos(4)-0.1;
set(gca,'OuterPosition', ca_pos);


hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend,'Location','northoutside','Orientation','horizontal','position',[0.5 0.95 0 0]);
legend('hide')

saveas(gcf, 'sim_experiment_planar_sigma' ,'epsc')
end
%%
function odinary_3d_results_time()
load ordinary3DresultsTime
figure('outerposition',[0    260    850    450])
yrange= [0 1100];

grid on
xdrawgraph(npts,yrange,method_list,'mean_c','Mean Cost',...
    'Number of Points','Cost (ms)');
hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend,'Location','northwest');
saveas(gcf, 'sim_experiment_ordinary_time' ,'epsc')
end
