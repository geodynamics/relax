% invert GPS data for static slip on multiple fault segments

clear all
close all
%clc

%% Data loading

settings=struct(...
    'Save_figures',false ...
    );

% Names of some of the fault segments
filename = 'coseismic_lowercrust_mantlelith';
if settings.Save_figures
    if 0==exist(['./' filename])
        mkdir(['./' filename])
    end
end

% Fault segments used in inversion
segment={...
    struct('name',['F2_' filename]), ...
    struct('name',['F3_' filename]), ...
    struct('name',['F4_' filename]), ...
    struct('name','yuha_30x48_305'), ...
    };

% sample segments to patches
fm=[];     % fault models
for i=1:length(segment)
    segment{i}.flt=[];
    fname=['./segments/' segment{i}.name '_seg.flt'];
                [local_index,x1,x2,x3,len,width,strike,dip,rake,lo,wo,alphal,alphaw,drake]=...
                    textread(fname,'%u %f %f %f %f %f %f %f %f %f %f %f %f %f%*[^\n]',...
                    'commentstyle','shell');
        for k=1:length(local_index)
            flt=flt2flt([x1(k);x2(k);x3(k)],len(k),width(k),strike(k),dip(k),rake(k),lo(k),wo(k),alphal(k),alphaw(k));
            flt=[flt,drake(k)*ones(size(flt,1),1)];
            segment{i}.flt=[segment{i}.flt;flt];
        end
    fm=[fm;segment{i}.flt];
end

% GPS data
gps=struct('network',['../GPS/postseismic_fit_5yr.dat'],'staNam',[],...
    'x',[],'y',[],'vx',[],'vy',[],'vz',[],'vxt',[],'vyt',[],'vzt',[], ...
    'd',[],'sig',[],'r',[],'G0',[],'G',[], ...
    'dt',[],'theta',[],'VECSIZE',3,'weight',[1 1 1]);

fid=fopen(gps.network);
gpsdata = textscan(fid,'%s %f %f %f %f %f %f %f %f','commentstyle','#');
fclose(fid);

gps.x = gpsdata{2};
gps.y = gpsdata{3};
gps.vx = 1e3*gpsdata{4}; %in mm
gps.vy = 1e3*gpsdata{5}; %in mm
gps.vz = 1e3*gpsdata{6}; %in mm
gps.vxsigma = 1e3*gpsdata{7};
gps.vysigma = 1e3*gpsdata{8};
gps.vzsigma = 1e3*gpsdata{9};

gps.d=zeros(gps.VECSIZE*length(gpsdata{1}),1); %initialize master interlaced GPS displacement vector
gps.sig=zeros(gps.VECSIZE*length(gpsdata{1}),1); %initialize master interlaced GPS error vector
gps.d(1:gps.VECSIZE:end)=gps.vx*gps.weight(1);
gps.sig(1:gps.VECSIZE:end)=gps.vxsigma*gps.weight(1);
gps.d(2:gps.VECSIZE:end)=gps.vy*gps.weight(2);
gps.sig(2:gps.VECSIZE:end)=gps.vysigma*gps.weight(2);
gps.d(3:gps.VECSIZE:end)=gps.vz*gps.weight(3);
gps.sig(3:gps.VECSIZE:end)=gps.vzsigma*gps.weight(3);
gps.Wx = diag(1./gps.sig.^2); %diagonal variance matrix Wx

% build GPS Green's function

% Poisson's solid
nu=1/4;

DGF=2;
rslip = [+1,+1];

for j=1:DGF
    G0=zeros(length(gps.d),size(fm,1));
    for k=1:size(fm,1)
        xg=fm(k,2);
        yg=fm(k,1);
        zg=fm(k,3);
        L=fm(k,4);
        W=fm(k,5);
        strike=fm(k,6)/180*pi;
        dip=fm(k,7)/180*pi;
        
        xd=gps.x-xg;
        yd=gps.y-yg;
        
        % strike- and dip-slip (j loop)
        [uxs,uys,uzs]=calc_okada(rslip(j),xd,yd,nu,dip,zg,L,W,j,strike);
        u = [uxs;uys;uzs]';
        
        % interlace east and north forward models
        for l=1:gps.VECSIZE
            G0(l:gps.VECSIZE:end,k)=gps.weight(l)*u(:,l);
        end
    end
    gps.G=[gps.G G0];
end

%% Smoothing matrix

gps.S = diag(gps.G'*gps.Wx*gps.G); %sensitivity
Sinvsqrt = 1./sqrt(gps.S); %S^-(1/2)

L=[];
% resolution-based smoothing
for i=1:length(segment)
    [xp,yp,zp,~]=transform4patch_general(segment{i}.flt(:,2),segment{i}.flt(:,1),segment{i}.flt(:,3),zeros(1,size(segment{i}.flt,1)),...
    segment{i}.flt(:,4),segment{i}.flt(:,5),segment{i}.flt(:,7),segment{i}.flt(:,6));
    Ltemp=compute_laplacian(mean(xp)',mean(yp)',mean(zp)',4);
    L=[L,zeros(size(L,1),size(Ltemp,2));
        zeros(size(Ltemp,1),size(L,2)),Ltemp];
    clear Ltemp
end

maxS = 4.4961e-04; %highest value of sensitivity
SinvsqrtT = maxS*[diag(Sinvsqrt(1:size(fm,1)))*L,zeros(size(L));zeros(size(L)),diag(Sinvsqrt(size(fm,1)+1:end))*L];

%% inversion

%************************
%       RHS VECTOR h
%************************
d=[gps.d];

% constraints
h=[d;
    zeros(size(SinvsqrtT,1),1)];

%************************
%       LHS MATRIX H
%************************
H=[gps.G;
    SinvsqrtT,zeros(size(SinvsqrtT,1),size(gps.G,2)-size(SinvsqrtT,2))];

%************************
%   POSITIVITY MATRIX C
%************************
rake1=fm(:,8)-fm(:,9);
rake2=fm(:,8)+fm(:,9);
C=-[diag(-sind(rake1)),diag(cosd(rake1));
    diag(sind(rake2)),diag(-cosd(rake2))];

options = optimoptions(@lsqlin,'MaxIter',1e+12,'Algorithm','active-set','Display','iter');

mt=lsqlin(H,h,C,zeros(size(C,1),1),[],[],[],[],[],options); % with positivity

%% GPS forward model

gps.dz = zeros(length(gps.d),1);
gps.dz(3:gps.VECSIZE:end) = gps.d(3:gps.VECSIZE:end);

gps.dt=H(1:length(gps.d),:)*mt;
gps.dtz = zeros(length(gps.dt),1);
gps.dtz(3:gps.VECSIZE:end) = gps.dt(3:gps.VECSIZE:end);

gps.dh=(gps.d - gps.dz);
gps.dth=(gps.dt - gps.dtz);
% residuals
gps.r=(gps.d-gps.dt);
gps.rx = (gps.d(1:gps.VECSIZE:end)-gps.dt(1:gps.VECSIZE:end));
gps.ry = (gps.d(2:gps.VECSIZE:end)-gps.dt(2:gps.VECSIZE:end));
gps.rz = (gps.d(3:gps.VECSIZE:end)-gps.dt(3:gps.VECSIZE:end));
gps.rh = gps.dh - gps.dth;

gps.vxt=gps.dt(1:gps.VECSIZE:end)/gps.weight(1);
gps.vyt=gps.dt(2:gps.VECSIZE:end)/gps.weight(2);
gps.vzt=gps.dt(3:gps.VECSIZE:end)/gps.weight(3);

% variance reduction
gps.theta=1-gps.r'*gps.r/(gps.d'*gps.d);
fprintf('   GPS variance reduction: %3.2f%% for ''%s''\n',gps.theta*100,gps.network)
gps.thetah=1-gps.rh'*gps.rh/(gps.dh'*gps.dh);
fprintf('   GPS variance reduction in h: %3.2f%% for ''%s''\n',gps.thetah*100,gps.network)
gps.thetax=1-gps.rx'*gps.rx/(gps.d(1:gps.VECSIZE:end)'*gps.d(1:gps.VECSIZE:end));
fprintf('   GPS variance reduction in x: %3.2f%% for ''%s''\n',gps.thetax*100,gps.network)
gps.thetay=1-gps.ry'*gps.ry/(gps.d(2:gps.VECSIZE:end)'*gps.d(2:gps.VECSIZE:end));
fprintf('   GPS variance reduction in y: %3.2f%% for ''%s''\n',gps.thetay*100,gps.network)
gps.thetaz=1-gps.rz'*gps.rz/(gps.d(3:gps.VECSIZE:end)'*gps.d(3:gps.VECSIZE:end));
fprintf('   GPS variance reduction in z: %3.2f%% for ''%s''\n',gps.thetaz*100,gps.network)

%% map view of surface uplift and subsidence

sx=201;
dx=2;
x=repmat((-sx/2:sx/2-1)'*dx,1,sx);
y=repmat((-sx/2:sx/2-1) *dx,sx,1);

Gup=zeros(sx*sx,DGF*size(fm,1));
for j=1:DGF
    for k=1:size(fm,1)
        
        xg=fm(k,2);
        yg=fm(k,1);
        zg=fm(k,3);
        Length=fm(k,4);
        W=fm(k,5);
        strike=fm(k,6)/180*pi;
        dip=fm(k,7)/180*pi;
        rake=fm(k,8)/180*pi;
        
        xd=x-xg;
        yd=y-yg;
        
        % Green's function
        [~,~,uzs]=calc_okada(rslip(j),xd,yd,nu,dip,zg,Length,W,j,strike); %Which calc_okada to use for 3D?
        
        % rake direction
        Gup(:,k+(j-1)*size(fm,1))=uzs';
    end
end

mtup=Gup*mt(1:DGF*size(fm,1));

mtup=reshape(mtup,sx,sx);
%% Plot results

cxlim=[-500 500];
cylim=[-500 500];
xlim = [-187.4 93.7];
ylim = [-83.3 166.5];

load ../plotting/ca_faults_dim.dat
load ../plotting/ca_faults_km.dat
load ../plotting/ca_coast_dim.dat
load ../plotting/ca_coast_km.dat
load ../plotting/ca_border_dim.dat
load ../plotting/ca_border_km.dat

load ../plotting/uplift_colormap.mat

fdurango=fopen(['../Wei_slipmodel/elmayor_wei_km_highslip.flt']);
durango = textscan(fdurango,'%u %f %f %f %f %f %f %f %f %f','commentstyle','#');
fclose(fdurango);

dflt=[durango{3},durango{4},durango{5},durango{6},durango{7},durango{8},durango{9},durango{10}];
dslip = 1e3*durango{2}; %in mm

%Import fault trace
ftrace = fopen(['../plotting/durango_trace_km.xyz']);
trace = textscan(ftrace,'%f %f %f');
fclose(ftrace);
tracex = trace{1};
tracey = trace{2};

%% Plot GPS and INSAR data

fontsize = 8;

figure(1); clf;
hold on
box on

[xp,yp,zp,up]=transform4patch_general(fm(:,2),fm(:,1),fm(:,3),sqrt(mt(1:size(fm,1)).^2 + mt(size(fm,1)+1:2*size(fm,1)).^2),fm(:,4),fm(:,5),fm(:,7),fm(:,6));
xpm=mean(xp);ypm=mean(yp);zpm=mean(zp);
patch(xp,yp,0*zp,up);

[xp,yp,~,~]=transform4patch_general(dflt(:,2),dflt(:,1),dflt(:,3),dslip*0,...
    dflt(:,4),dflt(:,5),dflt(:,7),dflt(:,6));
patch(xp,yp,0*xp,'EdgeColor',[0.45 0.2 0],'FaceColor','none','linewidth',0.5);

plot_faults(ca_faults_km,ca_faults_dim,xlim,ylim);
plot_faults(ca_coast_km,ca_coast_dim,cxlim,cylim);
plot_faults(ca_border_km,ca_border_dim,cxlim,cylim);

plot(tracex,tracey,'k','LineWidth',2);

axis equal tight
colormap([[linspace(1,0.9,32)';linspace(0.9,0.2,64)'],[linspace(1,0.6,32)';linspace(0.6,0.1,64)'],[linspace(1,0.2,32)';linspace(0.2,0.05,64)']])
set(gca,'clim',[0 2000]);
set(gca,'xlim',[-187.4 93.7],'ylim',[-55.5 166.5])

colorbar
xlabel('east (km)')
ylabel('north (km)')

if settings.Save_figures
    saveas(gcf,['./' filename '/slip_bigmapview.pdf'])
end

%%

figure(2); clf;
hold on
box on

scatter(x(:),y(:),25,mtup(:),'square','filled')

xlabel('east (km)')
ylabel('north (km)')
axis equal tight; box on

colormap(scec_colormap)

set(gca,'clim',[-20 20]);
set(gca,'xlim',[-187.4 93.7],'ylim',[-55.5 166.5])
colorbar('YTickLabel',{'-2.0 cm','-1.5 cm','-1.0 cm','-0.5 cm','0.0 cm','+0.5 cm','+1.0 cm','+1.5 cm','+2.0 cm'})
xlabel('east (km)')
ylabel('north (km)')

title 'GPS and INSAR data for first 6 months after earthquake'

if settings.Save_figures
    saveas(gcf,['./' filename '/verticals_mapview.pdf'])
end

%%

figure(3); clf;
hold on
box on

sc=0.5;
scatter(gps.x,gps.y,75,gps.vz,'filled')

quiver(gps.x,gps.y,sc*gps.vx ,sc*gps.vy ,0,'k','linewidth',1)

quiver(gps.x,gps.y,sc*gps.vxt,sc*gps.vyt,0,'color',[0 0.5 0],'linewidth',1)

plot_faults(ca_faults_km,ca_faults_dim,xlim,ylim);
plot_faults(ca_coast_km,ca_coast_dim,cxlim,cylim);
plot_faults(ca_border_km,ca_border_dim,cxlim,cylim);

xlabel('east (km)')
ylabel('north (km)')
axis equal tight; box on

colormap(scec_colormap)
set(gca,'clim',[-20 20]);
colorbar('YTickLabel',{'-2.0 cm','-1.5 cm','-1.0 cm','-0.5 cm','0.0 cm','+0.5 cm','+1.0 cm','+1.5 cm','+2.0 cm'})

set(gca,'xlim',[-187.4 93.7],'ylim',[-55.5 166.5])
xlabel('east (km)')
ylabel('north (km)')

title 'GPS and INSAR data for first 6 months after earthquake'

if settings.Save_figures
    saveas(gcf,['./' filename '/horizontals_mapview.pdf'])
end

%% plot fault model - total slip
for j=1:2
    figure(3+j);clf;
    set(gca,'xlim',[-70.3 70.3],'ylim',[-83.3 55.55],'zlim',[-25 0])
    
    xlim = [-70.3 70.3];
    ylim = [-83.3 55.5];
    colormap(scec_colormap)
    
    hold on
    box on
    
    % slip model = mt ; resolution = Rd2 ; bais = Rc2 ; discarded models = penalty
    
    switch j
        case 1
            toplot=mt(1:size(fm,1));
            title('JOINT INVERSION: modeled right-lateral slip');
        case 2
            toplot=mt(size(fm,1)+1:2*size(fm,1));
            title('JOINT INVERSION: modeled normal slip');
    end
    
    [xp,yp,zp,~]=transform4patch_general(fm(:,2),fm(:,1),fm(:,3),-toplot,...
        fm(:,4),min(fm(:,5),35),fm(:,7),fm(:,6));
    xpm=mean(xp);ypm=mean(yp);zpm=mean(zp);
    
    [~,~,~,up]=transform4patch_general(fm(:,2),fm(:,1),fm(:,3),-toplot,...
        fm(:,4),min(fm(:,5),35),fm(:,7),fm(:,6));
    
    pos=find(zpm>0);
    patch(xp(:,pos),yp(:,pos),-zp(:,pos),up(:,pos));
    
    [xp,yp,zp,up]=transform4patch_general(dflt(:,2),dflt(:,1),dflt(:,3),dslip,...
        dflt(:,4),dflt(:,5),dflt(:,7),dflt(:,6));
    %patch(xp,yp,-zp,up,'FaceColor','None','linewidth',0.5,'EdgeColor',[0.45 0.2 0]);
    
    %colormap([[linspace(0.5,1,96)';linspace(1,0.9,32)';linspace(0.9,0.2,64)'],[linspace(0,1,96)';linspace(1,0.6,32)';linspace(0.6,0.1,64)'],[linspace(0.75,1,96)';linspace(1,0.2,32)';linspace(0.2,0.05,64)']])
    colormap([[linspace(1,0.9,32)';linspace(0.9,0.2,64)'],[linspace(1,0.6,32)';linspace(0.6,0.1,64)'],[linspace(1,0.2,32)';linspace(0.2,0.05,64)']])
    
    axis equal tight
    set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-50 0])
    switch j
        case 1
            set(gca,'clim',[0 2000])
        case 2
            set(gca,'clim',[0 2000])
    end
    set(gca,'view',[0,0])
    colorbar
    
    if settings.Save_figures
        saveas(gcf,['./' filename '/Figure' num2str(4+j) '_short_joint.pdf'])
    end
end

%break

%% plot resolution of model parameters
figure(6);clf;
cla;
hold on

[xp,yp,zp,up]=transform4patch_general(fm(:,2),fm(:,1),fm(:,3),log10(gps.S(1:size(fm,1))),...
    fm(:,4),min(fm(:,5),35),fm(:,7),fm(:,6));

patch(xp,yp,-zp,up);

axis equal tight
box on
%colormap(cmap)
colormap(rot90(hot)')
colorbar
axis equal tight;
set(gca,'xlim',xlim,'ylim',ylim,'zlim',[-50 0])
set(gca,'clim',[-7 -4])
set(gca,'view',[0,0])

%title(['SVD estimate of resolution; ' num2str(M) ' parameters, trace(R)=' sprintf('%2.1f',sum(Restd(1:M)))])
title('JOINT INVERSION: SVD estimate of resolution')
xlabel('east (km)')
ylabel('north (km)')

if settings.Save_figures
    saveas(gcf,['./' filename '/Figure7_short_joint.pdf'])
end

%%

figure(7); clf;
hold on
box on

[xp,yp,zp,up]=transform4patch_general(fm(:,2),fm(:,1),fm(:,3),sqrt(mt(1:size(fm,1)).^2 + mt(size(fm,1)+1:2*size(fm,1)).^2),fm(:,4),fm(:,5),fm(:,7),fm(:,6));
xpm=mean(xp);ypm=mean(yp);zpm=mean(zp);
patch(xp,yp,0*zp,up);

[xp,yp,~,~]=transform4patch_general(dflt(:,2),dflt(:,1),dflt(:,3),dslip*0,...
    dflt(:,4),dflt(:,5),dflt(:,7),dflt(:,6));
patch(xp,yp,0*xp,'EdgeColor',[0.45 0.2 0],'FaceColor','none','linewidth',0.5);

plot_faults(ca_faults_km,ca_faults_dim,xlim,ylim);
plot_faults(ca_coast_km,ca_coast_dim,cxlim,cylim);
plot_faults(ca_border_km,ca_border_dim,cxlim,cylim);

axis equal tight
colormap([[linspace(1,0.9,32)';linspace(0.9,0.2,64)'],[linspace(1,0.6,32)';linspace(0.6,0.1,64)'],[linspace(1,0.2,32)';linspace(0.2,0.05,64)']])
set(gca,'clim',[0 2000]);
set(gca,'xlim',[-70.3 70.3],'ylim',[-83.3 55.5])

colorbar
xlabel('east (km)')
ylabel('north (km)')

if settings.Save_figures
    saveas(gcf,['./' filename '/slip_smallmapview.pdf'])
end
