function x = solver(A,b)
% solver: Solves a system of linear equations
%   x = solver(A,b):
%   Solves a system of linear equations using the '\' operator
%
% input:
%   A   = Coefficient matrix
%   b   = RHS vector
% output: 
%   x   = Solution vector
%           
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022

    x = A\b(:);
end