# Project: Euler and Midpoint Methods for ODEs

## 📌 Overview
Implements and compares three numerical schemes in MATLAB for solving linear initial value problems:

- Explicit Euler
- Implicit Euler
- Implicit Midpoint

Test problem:  
y'(t) = f1(t) + f2(t) * y(t), y(0) = y0
with constants `f1(t) = 2.0`, `f2(t) = -1.0`, and initial condition `y0 = 5`.

---

## 🧮 Methods
- **Explicit Euler**  
  y_{n+1} = y_n + h [ f1(t_n) + f2(t_n) y_n ]

- **Implicit Euler**  
  y_{n+1} = ( y_n + h f1(t_{n+1}) ) / ( 1 - h f2(t_{n+1}) )

- **Implicit Midpoint**  
  y_{n+1} = ( y_n + h f1(t_{m}) + (h f2(t_m)/2) y_n ) / ( 1 - (h f2(t_m)/2) ),  
  where t_m = (t_n + t_{n+1}) / 2

---

## 📊 Results
- Explicit Euler: simple but less stable.  
- Implicit Euler: more stable, especially for stiff problems.  
- Implicit Midpoint: combines stability with improved accuracy.  

*(You can add plots of the three methods vs exact solution if available.)*

---

## 📂 Files
- `euler_methods.m` → MATLAB implementation (explicit, implicit, midpoint).  
- `README.md` → project description.  

---

## ✅ Takeaway
This project illustrates how different time-stepping schemes behave in terms of stability and accuracy when applied to simple linear ODEs.
