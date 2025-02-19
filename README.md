# Boundary Element Method (BEM) Solver

## Overview
This repository contains a **2D Boundary Element Method (BEM) solver** for solving the **Laplace equation** on a circular domain. The solver applies **Dirichlet boundary conditions** on one half of the boundary and **Neumann boundary conditions** on the other half.

The implementation follows the methodology described in the book by **Why-Teong Ang**.

## Features
- Uses **Boundary Integral Method** to solve the Laplace equation.
- Supports **Dirichlet and Neumann boundary conditions**.
- Computes **boundary element properties** (midpoints, normals, and lengths).
- Constructs and solves the system of equations using numerical integration.
- Evaluates the solution inside the domain.
- Visualizes the computed solution.

## Getting Started
### Prerequisites
Ensure you have **MATLAB** installed on your system. This implementation does not require additional toolboxes.

### Running the Solver
1. Clone the repository or download the files.
2. Open **MATLAB** and navigate to the project directory.
3. Run the `main.m` script:
   ```matlab
   main
   ```

## File Descriptions
### **Main Solver Files**
- **`main.m`** - Entry point of the solver. Calls functions to construct the boundary model, apply boundary conditions, solve the system, and compute interior values.
- **`solver.m`** - Solves the system of equations formulated by the BEM method.

### **Boundary Model Construction**
- **`bem_model.m`** - Initializes the boundary element model and calls supporting functions to compute element properties.
- **`get_elements.m`** - Generates boundary elements (currently supports circular domains).
- **`get_midpoints.m`** - Computes midpoints of the boundary elements.
- **`get_normals.m`** - Computes normal vectors of the boundary elements.

### **Boundary Conditions**
- **`apply_boundary_conditions.m`** - Assigns Dirichlet and Neumann conditions to boundary elements.
- **`assign_solution.m`** - Updates the solution values in the BEM structure.

### **System of Equations & Numerical Integration**
- **`construct_axb.m`** - Constructs the linear system of equations (A*x = B) using boundary integrals.
- **`IF1.m` & `IF2.m`** - Perform numerical integration over boundary elements.

### **Solution Computation**
- **`calculate_domain.m`** - Computes the solution inside the domain at evenly distributed points.
- **`sol_point.m`** - Evaluates the solution at specific points inside the domain.

### **Visualization**
- **`plot_solution.m`** - Plots the computed solution over the domain.

## Example Usage
```matlab
main
```
This will generate the solution and display the result.

## References
- Why-Teong Ang, "A Beginner's Course in Boundary Element Methods", Springer.


[![View bem on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://in.mathworks.com/matlabcentral/fileexchange/104565-bem)
