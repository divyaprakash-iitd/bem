function plot_solution(x,y,phi,bem)
% plot_solution: Plots the solution as a contour plot
% plot_solution(x,y,phi):
%   Plots the solution as a contour plot along with the boundary elements
%   and the element normals
% input: 
%   x   =  A matrix containing the x-coordinates of the points
%   y   =  A matrix containing the y-coordinates of the points
%   phi =  A matrix containing the values of phi
%   bem =  A structure representing a boundary element model
%          bem.nelem       = No. of elements
%          bem.boundary.x  = x-coordinates of left node of boundary elements
%          bem.boundary.y  = y-coordinates of left node of boundary elements
%          bem.mid.x       = x-coordinates of mid-point of boundary elements
%          bem.mid.y       = y-coordinates of mid-point of boundary elements
%          bem.lelem       = Length of the elements
%          bem.normal.x    = x-component of the normal vector of elements
%          bem.normal.y    = y-component of the normal vector of elements
%          bem.bc.phi      = Value of the Dirichlet BC at each element
%          bem.bc.dphi     = Value of the Neumann BC at each element
%          bem.bc.type     = Type of BC of ab element->1:Dirichlet,2:Neumann
%
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022

    hold on
    for i = 1:bem.nelem
        plot(bem.boundary.x(i),bem.boundary.y(i),'ro')
        plot([bem.boundary.x(i),bem.boundary.x(i+1)],[bem.boundary.y(i),bem.boundary.y(i+1)],'k--')
        plot(bem.mid.x(i),bem.mid.y(i),'rx')
        quiver(bem.mid.x(i),bem.mid.y(i),bem.normal.x(i),bem.normal.y(i),0.5,'b')
    end
    colormap(jet)
    N = size(x,1);
    contourf(x,y,phi,N,'edgecolor','none')
    colorbar
    title('\phi(x,y)')
    axis equal
end