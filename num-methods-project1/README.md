# Project 1 – Numerical Methods for Differential Equations

This project investigates finite difference methods for numerical differentiation and several schemes for solving linear ODEs, including Euler methods, implicit midpoint, and Runge–Kutta methods.  
All implementations are written in MATLAB, with a focus on error analysis, convergence rates, and stability comparisons.

---

## 📌 Project Components

### Numerical Differentiation (Finite Differences)
- **File:** `sciprog_a12.m`  
- Implemented the one-sided finite difference quotient D_h f(x0) = (f(x0+h) - f(x0)) / h.  
- Evaluated derivatives of:
  - f(x) = exp(x) at x0 = 0 (exact derivative = 1),  
  - f(x) = sin(x) at x0 = 0 (exact derivative = 1),  
  - f(x) = x^(3/2) at x0 = 0 (exact derivative = 0).  
- Computed the error e_h = |f'(x0) - D_h f(x0)| for step sizes h = 2^-n, n=0…10.  
- Produced log–log error plots and tabulated convergence orders.  
- Observed first-order convergence for smooth functions (exp, sin), while the non-smooth function x^(3/2) showed slower convergence.

---

### Implementation of Euler & Midpoint Schemes
- **File:** `sciprog_a13.m`  
- Implemented three one-step methods as MATLAB functions with signatures:
  - `y = explicitEuler(t, f1, f2, y0)`
  - `y = implicitEuler(t, f1, f2, y0)`
  - `y = implicitMidpoint(t, f1, f2, y0)`  
- Each function computes the discrete solution y(t_n) for linear scalar IVPs of the form y'(t) = f1(t) + f2(t)y(t).  
- Verified correctness on small test cases.

---

### Testing the Schemes on Scalar Problems
- **Files:** `sciprog_a14a.m`, `sciprog_a14b.m`  
- Applied explicit Euler, implicit Euler, and implicit midpoint to:
  1. y'(t) = -y(t), y(0)=1, exact solution y(t)=exp(-t).  
     - Computed numerical solutions for step sizes h = 1/L.  
     - Plotted exact solution and all three approximations in the same figure.  
     - Computed error in the maximum norm: e_h = max |y(t_n) - y_n|.  
     - Plotted error vs h on a log–log scale to observe convergence rates.  
  2. y'(t) = 2 - y(t), y(0)=5.  
     - Compared growth/decay behavior across the three methods.  
     - Highlighted how implicit methods better capture stability.  

---

### Extending to Systems of ODEs
- **File:** `sciprog_a15.m`  
- Generalized all three methods to vector problems y'(t) = f1(t) + f2(t) y(t),  
  where f1 is a vector and f2 is a matrix.  
- Implemented the update rules in matrix form and verified they reduce to the scalar case when d=1.

---

### Linear System Case Study
- **Files:** `sciprog_a16a.m`, `sciprog_a16b.m`  
- Applied the vector methods to the coupled system:
  - y1'(t) = -y2(t),  
  - y2'(t) = y1(t),  
  with initial values y1(0)=1, y2(0)=0.  
- Exact solution: y(t) = (cos(t), sin(t)).  
- Computed numerical solutions with different step sizes.  
- Produced:
  - Error plots e_h vs h for explicit/implicit Euler and midpoint.  
  - Phase plots (y1 vs y2) showing how explicit Euler spirals outward (unstable), while implicit methods stay closer to the unit circle.  

---

### Runge–Kutta Methods
- **Files:** `sciprog_a17.m`, `sciprog_a18.m`, `sciprog_a19.m`  
- Implemented a general explicit Runge–Kutta solver given Butcher tableau (A, b, c).  
- Extended to systems (vector-valued problems).  
- Tested with the classical 4th-order Runge–Kutta method (RK4).  
- Compared convergence rates of Euler vs RK4 on y'(t)=-y(t).  
- Produced log–log error plots confirming 4th-order convergence of RK4.  
- Added phase plots for the system y1'=-y2, y2'=y1, showing RK4 preserves the circular trajectory much better than Euler methods.

---

## 📊 Key Outcomes
- Confirmed first-order convergence for the one-sided finite difference scheme on smooth functions.  
- Demonstrated through plots how explicit Euler is unstable for oscillatory systems, while implicit methods provide stability.  
- Showed midpoint and RK4 methods balance stability and accuracy, with RK4 achieving clear 4th-order error decay.  
- Phase plots provided strong visual evidence of stability differences across methods.  

---

## 📂 File Index
- `sciprog_a12.m` → Numerical Differentiation (Finite Differences)  
- `sciprog_a13.m` → Implementations of Explicit/Implicit Euler & Midpoint  
- `sciprog_a14a.m`, `sciprog_a14b.m` → Tests on scalar ODEs  
- `sciprog_a15.m` → Extension to systems of ODEs  
- `sciprog_a16a.m`, `sciprog_a16b.m` → Case study: coupled oscillatory system  
- `sciprog_a17.m` → General Runge–Kutta implementation  
- `sciprog_a18.m` → Runge–Kutta for systems  
- `sciprog_a19.m` → RK4 convergence and phase plots  
