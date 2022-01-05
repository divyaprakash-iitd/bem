clear; clc; close all;

% Description: Solves the 2D Laplace equation on an elliptical surface
% using the Boundary Integral Method

%% Number of elements
nelem = 20;

%% Generate the boundary element model
bem = bem_model(nelem,'circle');

%% Apply boundary conditions
bem = apply_boundary_conditions(bem);

%% Construct the system of equations
[A, B] = construct_axb(bem);

%% Solve for the unknowns and store them
sol = solver(A,B);

%% Assign solution to boundary
bem = assign_solution(bem,sol);

%% Solve for the interior domain
[x,y,phi] = calculate_domain(bem,nelem*10);

%% Plot the solution
plot_solution(x,y,phi)

%% Functions

function plot_solution(x,y,phi)
    figure(1)
    colormap(jet)
    N = size(x,1);
    contourf(x,y,phi,N,'edgecolor','none')
    colorbar
    title('\phi(x,y)')
    axis equal
end

function [x,y,phi] = calculate_domain(bem,N)
% Calculates the unknowns at N points distributed uniformly inside the
% domain
    r = linspace(0,0.9,N);
    theta = linspace(0,2*pi,N);
    [R, THETA] = meshgrid(r,theta);
    x = R.*cos(THETA);
    y = R.*sin(THETA);
    phi = zeros(N,N);
    % Calculate the solution
    for iN = 1:N
        for jN = 1:N
            xi  = x(iN,jN);
            eta = y(iN,jN);
            phi(iN,jN) = sol_point(bem,xi,eta);
        end
    end
end

function sum = sol_point(bem,xi,eta)
% Calculates the solution at a single point (xi,eta) inside the domain
    sum = 0;
    for i=1:bem.nelem
        xk = bem.boundary.x(i); 
        yk = bem.boundary.y(i);
        nkx = bem.normal.x(i); 
        nky = bem.normal.y(i);
        L = bem.lelem(i);
        F1 = IF1(xi,eta,xk,yk,nkx,nky,L);
        F2 = IF2(xi,eta,xk,yk,nkx,nky,L);
        ss = bem.bc.phi(i)*F2-bem.bc.dphi(i)*F1;
        sum = sum + ss;
    end
end
% %% Plot the domain
% hold on
% for i = 1:nelem
%     plot(xb(i),yb(i),'ro')
%     plot([xb(i),xb(i+1)],[yb(i),yb(i+1)],'k--')
%     plot(xm(i),ym(i),'rx')
%     quiver(xm(i),ym(i),nx(i),ny(i),0.5,'b')
% end
% axis equal
function x = solver(A,b)
    %x = A\b(:);
    x = bicg(A,b(:));
end

function bem = assign_solution(bem,sol)
    bem.bc.phi((bem.bc.type == 1)) = sol(bem.bc.type == 1);
    bem.bc.dphi((bem.bc.type == 0)) = sol(bem.bc.type == 0);
end
function bem = apply_boundary_conditions(bem)
    
    % Initialize arrays
    for i = 1:bem.nelem
        if (i <= bem.nelem/2)
            bem.bc.type(i) = 0;
            bem.bc.phi(i)    = 1;
        else 
            bem.bc.type(i) = 0;
            bem.bc.phi(i)    = 2;
        end
    end
end

function bem = bem_model(nelem,shape)
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

function [xb, yb, lelem] = get_elements(shape,nelem)
% get_midpoints: Calculates the mid-points
% [xm, ym] = get_midpoints(x,y):
%   Calculate the mid-points of the given array
%       X------------o------------X
%   (x_k,y_k)    (x_m,y_m)   (x_k+1,y_k+1)
%
% input: 
%   x     = A vector of size (N+1) containing the x-coordinates of the left
%           node of the N elements
%   y     = A vector of size (N+1) containing the y-coordinates of the left
%           node of the N elements
% output:
%   xm     = A vector of size N containing the x-coordinates of the 
%            mid-point of the N elements
%   xm     = A vector of size N containing the x-coordinates of the 
%            mid-point of the N elements
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

function [xm, ym] = get_midpoints(x,y)
% get_midpoints: Calculates the mid-points
% [xm, ym] = get_midpoints(x,y):
%   Calculate the mid-points of the given array
%       X------------o------------X
%   (x_k,y_k)    (x_m,y_m)   (x_k+1,y_k+1)
%
% input: 
%   x     = A vector of size (N+1) containing the x-coordinates of the left
%           node of the N elements
%   y     = A vector of size (N+1) containing the y-coordinates of the left
%           node of the N elements
% output:
%   xm     = A vector of size N containing the x-coordinates of the 
%            mid-point of the N elements
%   xm     = A vector of size N containing the x-coordinates of the 
%            mid-point of the N elements

    xm = 0.5*(x(1:end-1) + x(2:end));
    ym = 0.5*(y(1:end-1) + y(2:end));
end

function [nx, ny] = get_normals(x,y)
% get_midpoints: Calculates normal vector (pointing outward) components of 
% the elements
% [nx, ny] = get_normals(x,y):
%   Calculate the normal components of the given array
%                    ^
%                    |
%                    |
%                    |
%       X------------o------------X
%   (x_k,y_k)               (x_k+1,y_k+1)
%
% input: 
%   x     = A vector of size (N+1) containing the x-coordinates of the left
%           node of the N elements
%   y     = A vector of size (N+1) containing the y-coordinates of the left
%           node of the N elements
% output:
%   nx     = A vector of size N containing the x-components of the normal
%            vector of the N elements
%   ny     = A vector of size N containing the y-components of the normal
%            vector of the N elements

    l   = vecnorm([(x(2:end)-x(1:end-1)); (y(2:end)-y(1:end-1))]);
    ny  = (x(1:end-1)-x(2:end))./l;
    nx  = (y(2:end)-y(1:end-1))./l;   
end

function f = IF1(xi,eta,xk,yk,nkx,nky,L)
% IF1: Performs integration over an element
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

    A = L^2;
    B = (-nky*(xk-xi) + (yk-eta)*nkx)*2*L;
    E = (xk-xi)^2 + (yk-eta)^2;
   
    fun = @(t) log(A*t.^2 + B*t + E);
    
    f = L/4/pi*integral(fun,0,1);
end

function f = IF2(xi,eta,xk,yk,nkx,nky,L)
% IF1: Performs integration over an element
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
    A = L^2;
    B = (-nky*(xk-xi) + (yk-eta)*nkx)*2*L;
    E = (xk-xi)^2 + (yk-eta)^2;
    fun = @(t) (nkx*(xk-xi)+nky*(yk-eta))./(A*t.^2 + B*t + E);
    
    f = L/2/pi*integral(fun,0,1);
end