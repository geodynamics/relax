function [flt,wn,lt]=flt2flt(xo,L,W,strike,dip,rake,lo,wo,alphal,alphaw)
% [flt]=flt2flt(xo,L,W,strike,dip,lo,wo,alphal,alphaw)
%
% function flt2flt subsamples a fault patch into smaller segments
% of length and width starting from lo and wo, respectively
% increasing geometrically with down-dip distance with increment 
% alphal and alphaw (alphal>1 for increase).
%
% input:
%   xo     origin position vector [x1 (north);x2 (east);x3 (down)]
%   L      total length
%   W      total width
%   strike strike angle in degrees
%   dip    dip angle in degrees
%   rake   rake angle of slip
%   lo     approximative initial length of output segments
%   wo     approximative initial width of output segments
%   alpha1 geometric factor for length increase
%   alphaw geometric factor for width increase
%
% output:
%   flt    list of output segments in the format
%          x1,x2,x3,length,width,strike,dip
%
% to write an ascii output to use with the relax series:
%
%   fprintf('%3i %f %f %f %f %f %f %f %i\n',[[1:length(flt)]',flt]');

% create wi
Wc=W;
k=0;
w=0;
while Wc>0
    Wt=wo*alphaw^k;
    if Wt > Wc/2
        Wt = Wc;
    end
    wn=min([Wt,Wc]);
    w=[w; wn];
    k=k+1;
    Wc=Wc-wn;
end
Nw=k;

% strike and dip direction normal vectors
Sv=[ cosd(strike); sind(strike); 0];
Dv=[-cosd(dip)*sind(strike); cosd(dip)*cosd(strike); sind(dip)];

flt=[];

% loop in dip direction
k=0;
for j=1:Nw
    lt=lo*alphal^(j-1);
    Nl=ceil(L/lt);
    lt=L/Nl;
    
    % loop in strike direction
    for i=1:Nl
        x=xo+(i-1)*lt*Sv+sum(w(1:j))*Dv;
        k=k+1;
        slip=rand(1,1);
        flt=[flt; [x',lt,w(j+1),strike,dip,rake]]; 
    end
end