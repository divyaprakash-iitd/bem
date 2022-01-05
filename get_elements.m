function [xb, yb, lelem] = get_elements(shape,nelem)
%get_elements: Calculates the coordinates of the elements of a given shape
% [xb, yb, lelem] = get_elements(shape,nelem):
%   Calculates the coordinates of the elements of a given shape.
%   (Currently only 'circle' is supported)
%
% input: 
%   shape = A string specifiying the shape of the boundary.
%           (Currently only 'circle' is supported)
%   nelem = No. of boundary elements
% output:
%   xb     = A vector of size nelem+1 containing the x-coordinates of the 
%            left node of the nelem elements
%   yb     = A vector of size nelem+1 containing the y-coordinates of the 
%            left node of the nelem elements
%
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 05 January 2022

   xb = zeros(1,nelem+1);
   yb = zeros(1,nelem+1);
   
   if strcmp(shape,'circle')
       r = 1;
       theta = 0;
       dtheta = 2*pi/nelem;
       for i = 1:nelem
           xb(i) = r*cos(theta);
           yb(i) = r*sin(theta);
           theta = theta + dtheta;
       end
       
       xb(end) = xb(1);
       yb(end) = yb(1);
       
       lelem = vecnorm([(xb(2:end)-xb(1:end-1)); (yb(2:end)-yb(1:end-1))]);
   end
end
