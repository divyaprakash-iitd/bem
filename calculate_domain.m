function [x,y,phi] = calculate_domain(bem,N)
% calculate_domain: Calculates the unknowns at N^2 points distributed 
% uniformly inside the domain
%   [x,y,phi] = calculate_domain(bem,N):
%   Calculates the unknowns at N^2 points distributed uniformly inside the
%   domain
%
% input:
%   bem  =  A structure representing a boundary element model
%           bem.nelem       = No. of elements
%           bem.boundary.x  = x-coordinates of left node of boundary elements
%           bem.boundary.y  = y-coordinates of left node of boundary elements
%           bem.mid.x       = x-coordinates of mid-point of boundary elements
%           bem.mid.y       = y-coordinates of mid-point of boundary elements
%           bem.lelem       = Length of the elements
%           bem.normal.x    = x-component of the normal vector of elements
%           bem.normal.y    = y-component of the normal vector of elements
%           bem.bc.phi      = Value of the Dirichlet BC at each element
%           bem.bc.dphi     = Value of the Neumann BC at each element
%           bem.bc.type     = Type of BC of ab element->1:Dirichlet,2:Neumann
%   N   =   Square root of the number of points at which the solution
%           needs to be calculated
% output: 
%   x   =   A matrix containing the x-coordinates of the points
%   y   =   A matrix containing the y-coordinates of the points
%   phi =   A matrix containing the values of phi at the N^2 points
%           
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022

    r           = linspace(0,0.9,N);
    theta       = linspace(0,2*pi,N);
    [R, THETA]  = meshgrid(r,theta);
    x           = R.*cos(THETA);
    y           = R.*sin(THETA);
    phi         = zeros(N,N);
    
    % Calculate the solution
    for iN = 1:N
        for jN = 1:N
            xi  = x(iN,jN);
            eta = y(iN,jN);
            phi(iN,jN) = sol_point(bem,xi,eta);
        end
    end
end
