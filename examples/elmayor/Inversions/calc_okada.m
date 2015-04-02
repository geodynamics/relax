function [ux,uy,uz]=calc_okada(U,x,y,nu,delta,d,len,W,fault_type,strike);
%[ux,uy,uz]=calc_okada(U,x,y,nu,delta,d,len,W,fault_type,strike);
% U is slip
% x,y are the observation points
% d is depth of fault top (positive down)
% nu is Poisson ratio
% delta is dip angle (radians)
% len,W are the fault length and width, resp.
% fault_type is 1 2 3 for strike, dip, and opening
% uz is positive upwards
%
% the input is compatible with the convention of Wang and 
% Barbot (the RELAX software series).

x=reshape(x,1,numel(x));
y=reshape(y,1,numel(y));

cosd  = cos(delta);
sind  = sin(delta);

% define parameters with respect to the CENTER & BOTTOM of the fault (assume
% that the input ones correspond to the UPPER-LEFT TOP)

%s=[sin(strike) cos(strike) 0];
%d=[cos(strike)*cosd -sin(strike)*cosd sind];

d=d+W*sind;
x=x-W*cosd*cos(strike)-sin(strike)*len/2;
y=y+W*cosd*sin(strike)-cos(strike)*len/2;

strike = -strike+pi/2;
coss  = cos(strike);
sins  = sin(strike);
rot = [coss -sins ; sins coss];
rotx =  x*coss+y*sins;
roty = -x*sins+y*coss;

%%%%% Okada fault model for dislocation in an elastic half-space.
%%%%% based on BSSA Vol. 95 p.1135-45, 1985

L     = len/2;
Const = -U/(2*pi);
p = roty*cosd + d*sind;	%a matrix eqn. (30)
q = roty*sind - d*cosd;	%a matrix eqn. (30)
a = 1-2*nu;		% mu/(lambda+mu) = 1-2*poisson's ratio

parvec = [a, delta, fault_type];

[f1a,f2a,f3a] = fBi(rotx+L, p  , parvec, p, q);
[f1b,f2b,f3b] = fBi(rotx+L, p-W, parvec, p, q);
[f1c,f2c,f3c] = fBi(rotx-L, p  , parvec, p, q);
[f1d,f2d,f3d] = fBi(rotx-L, p-W, parvec, p, q);

%%%%% Displacement eqns. (25-27)

uxj = Const * (f1a - f1b - f1c + f1d);
uyj = Const * (f2a - f2b - f2c + f2d);
uz = Const * (f3a - f3b - f3c + f3d);

% rotate horizontals back to the orig. coordinate system
ux=-uyj*sins+uxj*coss;  
uy=uxj*sins+uyj*coss;

