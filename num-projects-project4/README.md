# Project 4 – Eigenvalue Iterations, 2D Quadrature, and FEM Data Structures

This project combines three major components:  
1. Iterative methods for eigenvalue computations.  
2. 2D quadrature rules on rectangles and triangles.  
3. Finite Element Method (FEM) data structures and efficient matrix assembly.  

All implementations are in MATLAB, using a mix of live scripts (`.mlx`) for demonstrations and functions (`.m`) for reusable components.

---

## 📌 Project Components

### Eigenvalue Iterations
- **Files:** `a41.mlx`, `a42.mlx`  
- Implemented and tested:
  - Power iteration → approximates the largest eigenvalue.  
  - Inverse iteration → approximates the smallest eigenvalue.  
  - Rayleigh quotient iteration → converges rapidly to nearby eigenvalues.  
- Compared against MATLAB’s `eigs`.  
- Analyzed sensitivity to the initial guess.  

---

### Condition Number Estimation
- **File:** `a42.mlx`  
- Used eigenvalue iterations to approximate the condition number κ₂(A) = λ_max / λ_min.  
- Tested on `gallery('poisson',k)` matrices with increasing size.  
- Verified scaling behavior with k².  

---

### 2D Quadrature on Rectangles
- **Files:** `a43.mlx`, `a44.mlx`  
- Extended 1D quadrature rules to 2D via tensor products.  
- Validated by integrating monomials x^j y^k exactly.  
- Applied Gauss quadrature to f(x,y)=exp(−(x²+y²)) over [0,1]×[0,2].  
- Observed exponential error decay in semilog plots.  

---

### 2D Quadrature on Triangles
- **Files:** `a45.mlx`, `a46.mlx`, `a47.mlx`, `a48.mlx`  
- Implemented quadrature rules on:
  - Reference triangle (Duffy transform).  
  - Arbitrary triangles via affine mappings.  
- Validated using monomials and smooth test functions f(x,y)=sin(2π(x+y)).  
- Showed exponential convergence until dominated by round-off errors.  

---

### FEM Data Structures (L-shape Domain)
- **Files:** `a49.mlx`, `Lshape_coordinates.dat`, `Lshape_elements.dat`, `Lshape_dirichlet.dat`  
- Parsed triangulation data for the L-shaped domain.  
- Visualized triangulation with interior and boundary edges.  
- Computed:
  - Domain area |Ω|.  
  - Boundary length.  

---

### FEM Quadrature on L-shape
- **File:** `a410.mlx`  
- Used triangle quadrature to approximate ∫_Ω f(x,y) dxdy for f(x,y)=exp(−(x²+y²)).  
- Repeated on refined meshes (`myrefine.m`).  
- Plotted quadrature error vs number of triangles N, showing algebraic convergence.  

---

### FEM Assembly and Efficiency
- **Files:** `a411.mlx`, `assemblyLaplace.m`, `assemblyLaplaceFast.m`, `provideGeometricData.m`, `myrefine.m`  
- Implemented assembly of FEM stiffness matrices:
  - Standard assembly (`assemblyLaplace.m`).  
  - Optimized sparse-aware assembly (`assemblyLaplaceFast.m`).  
- Benchmarked runtime for successive mesh refinements.  
- Log–log plots showed:
  - Naive assembly grows superlinearly.  
  - Sparse assembly scales nearly linearly.  

---

## 📊 Key Outcomes
- Demonstrated convergence and stability properties of eigenvalue iteration methods.  
- Verified condition number scaling for Poisson matrices.  
- Implemented and validated quadrature rules on 2D rectangles and triangles.  
- Built FEM data structures for the L-shape domain and computed geometric properties.  
- Showed convergence of quadrature-based FEM integrals with mesh refinement.  
- Benchmarked matrix assembly, confirming the efficiency of sparse data structures.  

---

## 📂 File Index
- `a41.mlx` → Eigenvalue iterations (power, inverse, Rayleigh)  
- `a42.mlx` → Condition number via eigenvalue iteration  
- `a43.mlx`, `a44.mlx` → 2D quadrature on rectangles  
- `a45.mlx`, `a46.mlx`, `a47.mlx`, `a48.mlx` → 2D quadrature on reference and arbitrary triangles  
- `a49.mlx` → FEM mesh data structure (L-shape domain)  
- `a410.mlx` → FEM quadrature on L-shape domain  
- `a411.mlx` → FEM assembly performance (naive vs sparse)  
- `assemblyLaplace.m`, `assemblyLaplaceFast.m` → FEM stiffness matrix assembly functions  
- `provideGeometricData.m` → FEM mesh preprocessing  
- `myrefine.m` → Mesh refinement utility  
- `Lshape_coordinates.dat`, `Lshape_elements.dat`, `Lshape_dirichlet.dat` → Triangulation data for L-shape  
