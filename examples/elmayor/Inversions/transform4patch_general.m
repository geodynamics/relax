function [xp,yp,zp,up]=transform4patch_general(x0,y0,z0,u0,dl,dw,dip,str)
% H. Perfettini
% S. Barbot : x0,y0,z0 is assumed top left corner
ntot=numel(x0);
if size(x0,2)==1,
    x0=x0';y0=y0';z0=z0';u0=u0';
    str=str';dip=dip';
end

x1=zeros(1,ntot);y1=zeros(1,ntot);z1=zeros(1,ntot);
x2=x1;y2=y1;z2=z1;
x3=x1;y3=y1;z3=z1;
x4=x1;y4=y1;z4=z1;

for kk=1:ntot,
    str0=str(kk);dip0=dip(kk);
    
    s=[sind(str0) cosd(str0) 0];
    d=[cosd(str0)*cosd(dip0) -sind(str0)*cosd(dip0) sind(dip0)];
    
    %building top left corner
    x1(kk)=x0(kk);
    y1(kk)=y0(kk);
    z1(kk)=z0(kk);
    %building top right corner
    x2(kk)=x0(kk)+dl(kk)*s(1);
    y2(kk)=y0(kk)+dl(kk)*s(2);
    z2(kk)=z0(kk);
    %building bottom right corner
    x3(kk)=x2(kk)+dw(kk)*d(1);
    y3(kk)=y2(kk)+dw(kk)*d(2);
    z3(kk)=z2(kk)+dw(kk)*d(3);
    %building bottom left corner
    x4(kk)=x1(kk)+dw(kk)*d(1);
    y4(kk)=y1(kk)+dw(kk)*d(2);
    z4(kk)=z1(kk)+dw(kk)*d(3);
end
xp=[x1(:) x2(:) x3(:) x4(:)]';
yp=[y1(:) y2(:) y3(:) y4(:)]';
zp=[z1(:) z2(:) z3(:) z4(:)]';
up=[u0(:) u0(:) u0(:) u0(:)]';