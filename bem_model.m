function bem = bem_model(nelem,shape)
% bem_model: Constructs a structure representing a boundary element model
% bem = bem_model(nelem,shape)
%   Constructs a structure representing a boundary element model for the
%   given shape
%
% input:
%   nelem = No. of boundary elements
%   shape = A string specifiying the shape of the boundary.
%           (Currently only 'circle' is supported)
%
% output: 
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
%           
%
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022

    [xb, yb, lelem]   = get_elements(shape,nelem);
    [xm, ym]          = get_midpoints(xb,yb);
    [nx, ny]          = get_normals(xb,yb);
    bem.nelem         = nelem;
    bem.boundary.x    = xb;
    bem.boundary.y    = yb;
    bem.mid.x         = xm;
    bem.mid.y         = ym;
    bem.lelem         = lelem;
    bem.normal.x      = nx;
    bem.normal.y      = ny;
    bem.bc.phi        = zeros(1,nelem);
    bem.bc.dphi       = zeros(1,nelem);
    bem.bc.type       = zeros(1,nelem);
end
