function [A, B] = construct_axb(bem)
% construct_axb: Constructs a system of linear equations of the form Ax=B
% [A, B] = construct_axb(bem):
%   Constructs a system of linear equations given a boundary element model
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
%           
% output:
%   A     = Coefficient matrix
%   B     = RHS vector
%
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022
    
    A = zeros(bem.nelem);
    B = zeros(1,bem.nelem);

    BCV = bem.bc.phi + bem.bc.dphi; 
    for m = 1:bem.nelem
        for k = 1:bem.nelem
            xk  = bem.boundary.x(k); 
            yk  = bem.boundary.y(k);
            xi  = bem.mid.x(m); 
            eta = bem.mid.y(m);
            nkx = bem.normal.x(k); 
            nky = bem.normal.y(k);
            L   = bem.lelem(k);
            F1  = IF1(xi,eta,xk,yk,nkx,nky,L);
            F2  = IF2(xi,eta,xk,yk,nkx,nky,L);

            delta = (m==k);

            % Fill up the LHS and RHS
            if (bem.bc.type(k) == 0)
                A(m,k)  = -F1;
                B(m)    = B(m) + BCV(k)*(-F2 + 0.5*delta);
            else
                A(m,k)  = F2-1/2*delta;
                B(m)    = B(m) + BCV(k)*F1;
            end 
        end
    end
end
