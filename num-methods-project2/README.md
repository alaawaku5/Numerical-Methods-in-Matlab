# Project 2 – Numerical Methods for Nonlinear ODEs

This project extends the study of numerical methods for ODEs from linear problems (Project 1) to nonlinear problems.  
The main focus is on fixed-point iteration, Newton’s method, implicit solvers for nonlinear ODEs, and adaptive step-size control.  
All implementations are written in MATLAB.

---

## 📌 Project Components

### Fixed-Point Iteration
- **File:** `sciprog_b21.m`  
- Implemented a general fixed-point solver `[z, zvec] = fixedpoint(G, z0)` for contractions G: R^d → R^d.  
- Used absolute error tolerance τabs = 1e-12 and maximum iterations Jmax = 1000.  
- Returned the final iterate and the sequence of iterates in zvec.  
- Verified linear convergence in test cases.

---

### Fixed-Point Test Problem
- **File:** `sciprog_b22.m`  
- Applied fixed-point iteration to the scalar contraction  
  G(z) = exp(-α|z|),  with α ∈ (0,1), initial guess z0=100.  
- Observed convergence to the unique fixed point z⋆ satisfying G(z⋆)=z⋆.  
- Demonstrated the effect of contraction mapping properties on convergence.

---

### One-Step Methods for Nonlinear ODEs
- **File:** `sciprog_b23.m`  
- Implemented three solvers for nonlinear ODEs y' = f(t,y):  
  1. **Explicit Euler** → direct update rule.  
  2. **Implicit Euler** → solved fixed-point equation z = yj + h f(tj+1, z) using `fixedpoint`.  
     - Returned also the iteration counts per step.  
  3. **Implicit Midpoint** → solved nonlinear equation  
     z = yj + h f((tj+1+tj)/2, (yj+z)/2).  
- Each solver returns the discrete trajectory matrix y(:,ℓ) plus, for implicit methods, iteration counts.

---

### Logistic Growth Model
- **File:** `sciprog_b24.m`  
- Considered the nonlinear logistic ODE:  
  y'(t) = k y(t)(C - y(t)),  with k=1, C=1, y0=0.5.  
- Exact solution: y(t) = 1 / (1 + (y0^{-1}-1)exp(-t)).  
- Applied Explicit Euler, Implicit Euler, and Implicit Midpoint on t ∈ [0,5].  
- Produced:  
  - Plot of exact vs numerical solutions.  
  - Plot of number of fixed-point iterations for implicit solvers.  
  - Error plots e_h vs h in log–log scale, confirming order of convergence.

---

### Newton’s Method for Nonlinear Systems
- **File:** `sciprog_b25.m`  
- Implemented Newton solver `[z, zvec] = newton(F, DF, z0)` with stopping criterion ∥z_{k+1}-z_k∥ ≤ τabs.  
- Verified quadratic convergence on test functions.  

---

### Newton Test Problems
- **File:** `sciprog_b26.m`  
- Tested Newton’s method on scalar nonlinear equations:  
  1. F(z)=z^2+exp(z)-2 → root ≈ 0.5373.  
  2. F(z)=log(z)+z^2 → root ≈ 0.6529.  
  3. F(z)=(cos(2z))^2 - z^2 → root ≈ 0.5149.  
- Initial guesses chosen in [0.25,10].  
- Plotted convergence orders q_k across iterations.  
- Observed quadratic convergence in successful runs.

---

### Lotka–Volterra Predator-Prey System
- **File:** `sciprog_b27.m`  
- Considered nonlinear system:  
  y1' = y1(1-y2),  
  y2' = -y2(1-y1),  
  with y(0) = (1,0.5).  
- Solved with Explicit Euler, Implicit Euler, and Implicit Midpoint up to T=15.  
- Compared:  
  - Fixed-point vs Newton iteration for implicit schemes (files `implicitEulerNewton`, `implicitMidpointNewton`).  
  - Plotted number of iterations per time step.  
- Showed Newton requires fewer iterations for convergence compared to fixed-point.

---

### Adaptive Time-Stepping (h–h/2 Strategy)
- **File:** `sciprog_b28.m`  
- Implemented function `[y, t] = adaptiveTimeStepping(solveODE, t0, T, f, y0, h0, tau)`.  
- Used the h–(h/2) strategy:  
  - Compute y_h with step h, y_h/2 with two half steps.  
  - Compare ε = ∥y_h - y_h/2∥∞ with tolerance τ h.  
  - Adapt step size accordingly (halve, keep, or double).  
- Ensured efficient step-size control while meeting accuracy requirements.

---

### Adaptive vs Uniform Stepping on Linear Oscillator
- **File:** `sciprog_b29.m`  
- Revisited the harmonic oscillator system from Project 1:  
  y1' = -y2, y2' = y1, with y(0) = (1,0).  
- Compared Explicit Euler with:  
  - Uniform mesh refinement (h=2^-n).  
  - Adaptive refinement (tolerance τ=2^-n, initial step size h0=5).  
- Plotted error at final time e_N = |y(2π)-y_N| vs number of steps N in log–log scale.  
- Observed convergence rate e_N = O(N^-α) and demonstrated efficiency of adaptive stepping.

---

## 📊 Key Outcomes
- Verified linear convergence of fixed-point iteration and quadratic convergence of Newton’s method.  
- Applied implicit solvers with nonlinear iteration to logistic growth and predator-prey systems.  
- Demonstrated efficiency gains of Newton’s method over fixed-point in nonlinear ODE solvers.  
- Implemented adaptive time stepping and showed it achieves target accuracy with fewer steps compared to uniform grids.  

---

## 📂 File Index
- `sciprog_b21.m` → Fixed-point solver  
- `sciprog_b22.m` → Fixed-point test problem  
- `sciprog_b23.m` → Explicit/Implicit Euler and Midpoint solvers (nonlinear ODEs)  
- `sciprog_b24.m` → Logistic growth model study  
- `sciprog_b25.m` → Newton’s method implementation  
- `sciprog_b26.m` → Newton test problems with convergence plots  
- `sciprog_b27.m` → Lotka–Volterra predator-prey system (fixed-point vs Newton)  
- `sciprog_b28.m` → Adaptive time-stepping function (h–h/2 strategy)  
- `sciprog_b29.m` → Adaptive vs uniform stepping on harmonic oscillator  
