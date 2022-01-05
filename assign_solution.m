
function bem = assign_solution(bem,sol)
% assign_solution: Assigns the calculated values at the boundaries to their
% respective arrays in the boundary element model
% bem = bem_model(nelem,shape)
%   Constructs a structure representing a boundary element model for the
%   given shape
%
% input:
%   bem = A structure representing a boundary element model
%         bem.nelem       = No. of elements
%         bem.boundary.x  = x-coordinates of left node of boundary elements
%         bem.boundary.y  = y-coordinates of left node of boundary elements
%         bem.mid.x       = x-coordinates of mid-point of boundary elements
%         bem.mid.y       = y-coordinates of mid-point of boundary elements
%         bem.lelem       = Length of the elements
%         bem.normal.x    = x-component of the normal vector of elements
%         bem.normal.y    = y-component of the normal vector of elements
%         bem.bc.phi      = Value of the Dirichlet BC at each element
%         bem.bc.dphi     = Value of the Neumann BC at each element
%         bem.bc.type     = Type of BC of ab element->1:Dirichlet,2:Neumann
%   sol = The solution obtained by solving the system of linear equations
% output: 
%   bem = A structure representing a boundary element model 
%           
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022

    bem.bc.phi((bem.bc.type == 1))  = sol(bem.bc.type == 1);
    bem.bc.dphi((bem.bc.type == 0)) = sol(bem.bc.type == 0);
end
