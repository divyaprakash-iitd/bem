function bem = apply_boundary_conditions(bem,phi,dphi)
% apply_boundary_conditions: Applies boundary conditions to provided 
% boundary element model
% bem = apply_boundary_conditions(bem,phi,dphi):
%   Applies boundary conditions to the boundary of the provided boundary
%   element model. Dirichlet and Neumann boundary conditions are applied to
%   each half of the boundary.
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
%   phi  =  Value of  Dirichlet boundary conditions
%   dphi =  Value of  Dirichlet boundary conditions
%
% output: 
%   bem  = A structure representing a boundary element model
%           
%
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022

    for i = 1:bem.nelem
        if (i <= bem.nelem/2)
            bem.bc.type(i) = 0;
            bem.bc.phi(i)  = phi;
        else 
            bem.bc.type(i) = 1;
            bem.bc.dphi(i)  = dphi;
        end
    end
end
