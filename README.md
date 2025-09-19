# Numerical Methods in MATLAB

This repository contains a collection of projects developed as part of the course *Scientific Programming for Interdisciplinary Mathematics*.  
Each project explores a different aspect of **numerical analysis** and **scientific computing**, implemented in MATLAB with a focus on accuracy, convergence, and efficiency.  

The projects progress from basic finite difference approximations to advanced finite element methods (FEM) in 2D.  
All codes are self-contained and include plots, error analysis, and comparisons with exact solutions where available.  

---

## 📂 Projects Overview

### [Project 1 – Numerical Methods for Differential Equations](num-methods-project1)
- Finite difference approximations of derivatives.  
- Explicit Euler, Implicit Euler, and Implicit Midpoint schemes.  
- Extension to systems of ODEs and phase plots.  
- Runge–Kutta methods (including RK4) with convergence studies.  

### [Project 2 – Nonlinear ODE Solvers](num-methods-project2)
- Fixed-point iteration and Newton’s method for nonlinear equations.  
- Nonlinear ODE solvers with implicit methods (using iteration).  
- Logistic growth and Lotka–Volterra predator-prey systems.  
- Adaptive time-stepping strategies (h–h/2 method).  

### [Project 3 – Quadrature and 1D FEM](num-methods-project3)
- Numerical quadrature: equidistant vs Gauss–Legendre nodes.  
- Error analysis for smooth vs singular functions.  
- Implementation of P1 finite element method in 1D.  
- Diffusion, advection, and reaction terms with varying coefficients.  

### [Project 4 – Eigenvalue Iterations and 2D FEM Preliminaries](num-methods-project4)
- Power, inverse, and Rayleigh quotient iterations for eigenvalues.  
- Condition number estimation for large sparse matrices.  
- 2D quadrature on rectangles and triangles (Duffy transform).  
- FEM data structures for the L-shaped domain.  
- Sparse vs naive assembly of FEM stiffness matrices.  

---

## 🛠 Tools & Methods
- **Language:** MATLAB  
- **Numerical techniques:** finite differences, explicit/implicit Euler, midpoint methods, Runge–Kutta, fixed-point iteration, Newton’s method, quadrature rules, FEM assembly.  
- **Focus areas:** error analysis, stability, convergence rates, computational efficiency.  

---

## 🎯 Key Learning Outcomes
- Verified theoretical convergence rates through experiments.  
- Visualized stability differences between explicit and implicit methods.  
- Implemented quadrature rules and FEM solvers from scratch.  
- Benchmarked sparse matrix assembly, confirming efficiency gains.  

---

📌 This repository serves as a portfolio of numerical methods implementations — useful for reference, teaching, and as a foundation for more advanced research in numerical PDEs and FEM.
