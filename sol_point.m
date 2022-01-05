function sol = sol_point(bem,xi,eta)
% sol_point: Solves for the solution at a point
%   sum = sol_point(bem,xi,eta):
%   Solves for the solution at a point (xi, eta)

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
%   xi  =   x-ccordinate of the point
%   eta =   x-ccordinate of the point
% output: 
%   sol = Solution at the provided point
%           
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022

    sol = 0;
    for i=1:bem.nelem
        xk = bem.boundary.x(i); 
        yk = bem.boundary.y(i);
        nkx = bem.normal.x(i); 
        nky = bem.normal.y(i);
        L = bem.lelem(i);
        F1 = IF1(xi,eta,xk,yk,nkx,nky,L);
        F2 = IF2(xi,eta,xk,yk,nkx,nky,L);
        ss = bem.bc.phi(i)*F2-bem.bc.dphi(i)*F1;
        sol = sol + ss;
    end
end
