clear all
close all
clc

fdata = '../../Papers/Heat flow/USGS_CaliforniaHeatFlowData_mod.txt';
fid = fopen(fdata);
input = textscan(fid,'%f %f %f');
fclose(fid);
norths = (input{1} - 32.5)*111.1;
easts = (input{2} + 115.5)*93.7;
HF = input{3};

scatter(easts,norths,50,HF,'+','LineWidth',1)
axis equal
xlim([-187.4 93.7]);
ylim([-55.5 166.5]);

load hot2.mat
colormap(hot2)
%colormap([[linspace(1,1,32),linspace(1,0.5,32)]',[linspace(1,1,32),linspace(1,0,32)]',[linspace(1,0,32),linspace(0,0,32)]'])
colorbar
hold on
plot(23.5,-11.1,'^','Color','k')
set(gca,'clim',[0 150])
box on

saveas(gcf,'../../Papers/heatflowmap.pdf')