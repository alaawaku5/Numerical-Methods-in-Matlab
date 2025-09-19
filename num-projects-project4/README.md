# Project 3 – Quadrature and 1D Finite Element Method

This project investigates **numerical quadrature (integration)** methods and the implementation of the **1D P1 Finite Element Method (FEM)** for elliptic boundary value problems.  
All code is written in MATLAB. The project highlights how quadrature accuracy affects PDE solvers and culminates in a full FEM implementation.

---

## 📌 Project Components

### Quadrature Weights via Linear Systems
- **File:** `a31.m`  
- Implemented computation of quadrature weights for given nodes on [0,1] by solving a linear system Aw=b.  
- Verified exact integration of polynomials up to degree n−1.

---

### Quadrature with Equidistant Nodes
- **File:** `a32.m`  
- Tested equidistant nodes x_j = j/n for n = 1,…,20.  
- Integrated:
  - f(x)=exp(x) → observed exponential error decay (semilog plot).  
  - f(x)=x^(3/2) → observed algebraic error decay (log–log plot).  

---

### Quadrature on Arbitrary Intervals
- **File:** `a33.m`  
- Implemented `[x, w] = quadrature(a, b, n)` that rescales quadrature from [0,1] to arbitrary [a,b].  
- Verified using substitution φ(t)=a+t(b−a).

---

### Testing Arbitrary Interval Quadrature
- **File:** `a34.m`  
- Computed ∫_{−π}^{exp(1)} exp(−x) dx ≈ 23.0747.  
- Compared error vs number of quadrature points n=1,…,20.  
- Confirmed convergence rate in semilog plots.

---

### Gauss–Legendre Quadrature
- **Files:** `a35.m`, `gaussLegendre.m`  
- Adapted Gauss–Legendre quadrature from [−1,1] to arbitrary [a,b].  
- Implemented `[x, w] = quadratureGauss(a, b, n)`.  
- Verified that Gauss quadrature integrates polynomials up to degree 2n−1 exactly.

---

### Testing Gauss Quadrature
- **File:** `a36.m`  
- Repeated the test integral from `a34.m` with Gauss quadrature.  
- Added error curves for Gauss quadrature vs equidistant quadrature.  
- Observed superior accuracy of Gauss quadrature.

---

### Advanced Quadrature Test
- **Files:** `a37_part1.m`, `a37_part2.m`  
- Computed ∫_{−5}^{5} 1/(1+x^2) dx ≈ 2.7468 with both equidistant and Gauss quadrature.  
- Compared accuracy for n=1,…,50 nodes.  
- Observed faster convergence of Gauss quadrature, especially for smooth integrands.

---

### Implementing P1 Finite Element Method (Diffusion Only)
- **File:** `a38.m`, `p1FEM.m`  
- Implemented FEM solver for −(a u′)′ = f with u(0)=u(1)=0 on [0,1].  
- Used hat functions as basis and numerical quadrature for assembling stiffness matrix A and load vector L.  
- Verified with a(x)=1, f(x)=1 against exact solution u(x)=½(1−x)x.

---

### General Second-Order Elliptic Problems
- **File:** `a39_310.m`  
- Extended FEM to include:
  - Diffusion term (A matrix)  
  - Advection term (B matrix)  
  - Reaction term (C matrix)  
- Implemented `x = p1FEM(a, b, c, f, t)` handling general coefficients.  
- Assembled full system Sx=L with contributions from A, B, and C.

---

### Diffusion–Advection–Reaction Test Cases
- **Files:** `a311.m`, `a311_2.m`  
- Tested FEM for a(x)=1, c(x)=−1, f(x)=1, and b(x) ∈ {−10,−5,0,5,10}.  
- Plotted solutions for different advection constants b.  
- Demonstrated how increasing |b| skews the solution profile, illustrating the role of advection.

---

### Efficiency and Sparse Matrices
- **File:** `last.m` (if used)  
- Investigated computational time for FEM assembly with N=2^n elements.  
- Produced log–log plots of runtime vs N, showing polynomial growth.  
- Discussed how FEM matrices are tridiagonal/sparse, and how MATLAB’s `sparse` format can reduce storage and improve efficiency.

---

## 📊 Key Outcomes
- Verified exponential vs algebraic convergence of quadrature for smooth vs singular integrands.  
- Demonstrated superior accuracy of Gauss quadrature compared to equidistant-node quadrature.  
- Implemented the full P1 FEM for 1D elliptic PDEs, including diffusion, advection, and reaction terms.  
- Showed the effect of advection on solution shape.  
- Identified efficiency gains possible with sparse matrix assembly in FEM.

---

## 📂 File Index
- `a31.m` → Quadrature weights (linear system formulation)  
- `a32.m` → Equidistant quadrature nodes, error plots  
- `a33.m` → Quadrature on arbitrary intervals  
- `a34.m` → Test of quadrature on ∫ exp(−x) dx over [−π,exp(1)]  
- `a35.m`, `gaussLegendre.m` → Gauss–Legendre quadrature implementation  
- `a36.m` → Gauss quadrature test vs equidistant quadrature  
- `a37_part1.m`, `a37_part2.m` → Quadrature test on ∫ 1/(1+x^2) dx  
- `a38.m`, `p1FEM.m` → FEM implementation for diffusion-only problems  
- `a39_310.m` → Extended FEM for diffusion–advection–reaction problems  
- `a311.m`, `a311_2.m` → Test cases with varying advection parameter b  
- *(optional)* `last.m` → Runtime scaling and sparse matrix discussion  
