function Lap = compute_laplacian(x,y,z,N)
%COMPUTE_LAPLACIAN   Generate Laplacian with respect to a set of points
%   LAP = COMPUTE_LAPLACIAN(X,Y,Z,N) takes in the coordinates of a set of
%   points that are randomly scattered on an unknown surface and generates
%   a discrete approximation of the Laplacian using the nearest N points.
%   X, Y, and Z are column vectors of the same length, and N must be a
%   positive integer.
%
%   Example
%   x = repmat([0:5],6,1); x = x(:);
%   y = repmat([0:5]',1,6); y = y(:);
%   z = repmat([0:0.1:0.5],6,1); z = z(:);
%   N = 4;
%   Lap = compute_laplacian(x,y,z,N)
%   
%   See also compute_laplacian, get_fault_model, PCAIM_driver.

%   By Andrew Kositsky
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.0.0.0 $  $Date: 2009/12/03  $


% Define the N_nearest variable
N_nearest = N;

% Define the number of points on the surface
n_points = numel(x);

% Define the data matrix, where each column is a different coordinate
data = [x,y,z];

% Pre-allocate the distance matrix and the Laplacian matrix
distance_matrix = zeros(n_points);
Lap = zeros(n_points);

% Find the distance between every two points on the lower triangular
% portion of the distance_matrix
for ii = 1:n_points-1
    for jj = ii+1 : n_points
        distance_matrix(ii,jj) = dist_fcn(data(ii,:),data(jj,:));
    end
end

%%% Fill opposite entries of matrix, relying on the fact that dist_fcn(x,x) = 0
distance_matrix = distance_matrix + distance_matrix';

% Sort the distances, keep the index, and pull out the indexes of the
% N_nearest nearest points to each point on the fault
[dummy_var,I_dist] = sort(distance_matrix,2); clear('dummy_var');
NearestNPoints = (I_dist(:,2:1+N_nearest));

% Calculate the Laplacian individually for each point. See:
%   Geertjan Huiskamp. Difference formulas for the surface laplacian on a
%   triangulated surface. Journal of Computational Physics, 95(2):477 -
%   496, 1991. 
%   or the PCAIM manual for details.

for ii = 1:n_points
    % center the data at data(ii)
    translated_data = data(NearestNPoints(ii,:),:)-repmat(data(ii,:),N_nearest,1);
    % find the best fitting plane going through data(ii), namely that
    % spanned by the columns of v. And the distance of the points'
    % projection on this plane distances along these axes are u * s.
    [u,s,v] = svds(translated_data,2);
    
    % Convert to polar coordiantes to sort by angle
    [theta,r] = cart2pol(u(:,1)*s(1,1),u(:,2)*s(2,2));
    
    % Sort by angle
    [theta,theta_index] = sort(theta);
    
    % Reorder the points for NearestNPoints, r, and u
    NearestNPoints(ii,:) = NearestNPoints(ii,theta_index);
    r = r(theta_index);
    u = u(theta_index,:);
    
    % Calculate the differences in theta for use in the Laplacian
    %    approximation
    delta_theta = diff([theta;theta(1)]);
    
    % Assign thetas in the positive rotation direction
    theta_plus = delta_theta;
    % Assign thetas in the negative rotation direction
    theta_minus = [delta_theta(end);delta_theta(1:end-1)];
    % Assign the average distance to the neighboring points
    r_bar = mean(r);
    
    % Calculate Theta_tot as in (Huiskamp, 1991)
    Theta_tot = calc_Theta_tot(theta_minus,theta_plus);
    
    % From page 4 of my Non-Planar Fault Laplacian 0.2 notes, equation 12;
    % or PCAIM manual.
    for jj = 1:N_nearest
        Lap(ii,NearestNPoints(ii,jj)) = 4/r_bar * 1/Theta_tot * 1/r(jj) *...
            calc_Theta(theta_plus(jj),theta_minus(jj));
    end
    % Central point is equal to the negative sum of the neighboring weights
    Lap(ii,ii) = - sum(Lap(ii,NearestNPoints(ii,:)));
end



function dist = dist_fcn(x1,x2)

%DIST_FCN   Find the Eucliean norm between two input vectors
dist = norm(x1-x2);

function Theta_tot = calc_Theta_tot(theta_plus,theta_minus)
%CALC_THETA_TOT   Defined as in (Huiskamp, 1991)
Theta_tot = 0;
for i = 1 : numel(theta_plus)
    Theta_tot = Theta_tot + calc_Theta(theta_plus(i),theta_minus(i));
end

function Theta = calc_Theta(theta_plus,theta_minus)
%CALC_THETA   Defined as in (Huiskamp, 1991)
Theta = calc_half_Theta(theta_plus) + calc_half_Theta(theta_minus);

function partial_Theta = calc_half_Theta(theta_diff)
%CALC_HALF_THETA   Defined in (Huiskamp, 1991) for theta_plus/theta_minus
tol = 10^(-10);

if(abs(theta_diff) < tol)
    partial_Theta = 0;
else
    partial_Theta = (1-cos(theta_diff))/sin(theta_diff);
end


function iscolinear_flag = iscolinear(x1,x2,x3)
%ISCOLINEAR   check if three points are colinear up to a tolerance.
%   This function is here to avoid dividing by zero when we look at
%   sin(theta), where theta is the angle between x1, x2, and x3.
%
tol = 10^(-10);

% if the dot product between the two vectors is very close to 1,
if(dot((x1-x2)/norm(x1-x2),(x1-x3)/norm(x1-x3)) > 1-tol)
    %they are colinear
    iscolinear_flag = 1;
else
    %else they are not
    iscolinear_flag = 0;
end