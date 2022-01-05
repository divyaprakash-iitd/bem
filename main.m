clear; clc; close all;

%% Description:
% Solves the 2D Laplace equation on a circular surface using
% using the Boundary Integral Method. The boundary conditions are specified
% as Dirichlet on one half of the boundary and as Neumann condition on the
% other half.

%% Number of elements
nelem = 7;

%% Generate the boundary element model
bem = bem_model(nelem,'circle'); % Currently only circle is supported

%% Apply boundary conditions
phi_v   = 1;
dphi_v  = 2;
bem     = apply_boundary_conditions(bem,phi_v,dphi_v);

%% Construct the system of equations
[A, B] = construct_axb(bem);

%% Solve for the unknowns and store them
sol = solver(A,B);

%% Assign solution to boundary
bem = assign_solution(bem,sol);

%% Solve for the interior domain
[x,y,phi] = calculate_domain(bem,nelem+1);

%% Plot the solution
plot_solution(x,y,phi,bem)
