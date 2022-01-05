function [xm, ym] = get_midpoints(x,y)
% get_midpoints: Calculates the mid-points
% [xm, ym] = get_midpoints(x,y):
%   Calculate the mid-points of the given array
%       X------------o------------X
%   (x_k,y_k)    (x_m,y_m)   (x_k+1,y_k+1)
%
% input: 
%   x     = A vector of size (N+1) containing the x-coordinates of the left
%           node of the N elements
%   y     = A vector of size (N+1) containing the y-coordinates of the left
%           node of the N elements
% output:
%   xm     = A vector of size N containing the x-coordinates of the 
%            mid-point of the N elements
%   xm     = A vector of size N containing the x-coordinates of the 
%            mid-point of the N elements

    xm = 0.5*(x(1:end-1) + x(2:end));
    ym = 0.5*(y(1:end-1) + y(2:end));
end