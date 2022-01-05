function [nx, ny] = get_normals(x,y)
% get_midpoints: Calculates normal vector (pointing outward) components of 
% the elements
% [nx, ny] = get_normals(x,y):
%   Calculate the normal components of the given array
%                    ^
%                    |
%                    |
%                    |
%       X------------o------------X
%   (x_k,y_k)               (x_k+1,y_k+1)
%
% input: 
%   x     = A vector of size (N+1) containing the x-coordinates of the left
%           node of the N elements
%   y     = A vector of size (N+1) containing the y-coordinates of the left
%           node of the N elements
% output:
%   nx     = A vector of size N containing the x-components of the normal
%            vector of the N elements
%   ny     = A vector of size N containing the y-components of the normal
%            vector of the N elements
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022

    l   = vecnorm([(x(2:end)-x(1:end-1)); (y(2:end)-y(1:end-1))]);
    ny  = (x(1:end-1)-x(2:end))./l;
    nx  = (y(2:end)-y(1:end-1))./l;   
end
