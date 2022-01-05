function f = IF2(xi,eta,xk,yk,nkx,nky,L)
% IF2: Performs integration over an element
% the elements
% f = IF1(xi,eta,xk,yk,nkx,nky,L):
%   Performs integration over an element
% input:
%   xi      = x location on element
%   eta     = y location on element
%   xk      = x-coordinate of the left node of the boundary element
%   yk      = y-coordinate of the left node of the boundary element
%   nkx     = x-component of normal vector of the element
%   nky     = y-component of normal vector of the element
%   L       = Length of the element
% output:
%   f       = Integral value
%
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022

    A = L^2;
    B = (-nky*(xk-xi) + (yk-eta)*nkx)*2*L;
    E = (xk-xi)^2 + (yk-eta)^2;
    fun = @(t) (nkx*(xk-xi)+nky*(yk-eta))./(A*t.^2 + B*t + E);
    
    f = L/2/pi*integral(fun,0,1);
end
