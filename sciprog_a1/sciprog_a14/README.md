# Project: Euler and Midpoint Methods for ODEs (Exercise 4)

## 📌 Overview
This project implements and compares three numerical schemes in MATLAB for solving linear initial value problems (IVPs):

- **Explicit Euler**
- **Implicit Euler**
- **Implicit Midpoint**

We apply these methods to two test problems:

### Part (a)
\[
y'(t) = -y(t), \quad y(0) = 1, \quad \text{Exact solution: } y(t) = e^{-t}.
\]

### Part (b)
\[
y'(t) = f_1(t) + f_2(t) y(t), \quad y(0) = 5,
\]
with constants \( f_1(t) = 2.0 \), \( f_2(t) = -1.0 \).

---

## 🧮 Methods
- **Explicit Euler**  
  \[
  y_{n+1} = y_n + h [ f_1(t_n) + f_2(t_n) y_n ]
  \]

- **Implicit Euler**  
  \[
  y_{n+1} = \frac{y_n + h f_1(t_{n+1})}{1 - h f_2(t_{n+1})}
  \]

- **Implicit Midpoint**  
  \[
  y_{n+1} = \frac{y_n + h f_1(t_m) + \tfrac{h f_2(t_m)}{2} y_n}{1 - \tfrac{h f_2(t_m)}{2}},
  \quad t_m = \tfrac{t_n + t_{n+1}}{2}
  \]

---

## 📊 Results
- **Part (a)**:  
  - Explicit Euler: conditionally stable, less accurate.  
  - Implicit Euler: stable but slightly diffusive.  
  - Implicit Midpoint: stable and more accurate.  

- **Part (b)**:  
  - Similar behavior observed with forcing term included.  
  - Implicit schemes outperform explicit Euler in terms of stability.  

### Example Output (Part a)
![Euler Methods Comparison](results_part_a.png)

### Example Output (Part b)
![Euler Methods Comparison](results_part_b.png)

---

## 📂 Files
- `exercise4.m` → MATLAB code (both parts a and b, with exact solution for part a).  
- `results_part_a.png` → plot comparing exact and numerical solutions for part a.  
- `results_part_b.png` → plot comparing explicit, implicit, and midpoint for part b.  

---

## ✅ Takeaway
This project highlights how explicit and implicit methods differ in accuracy and stability for linear ODEs.  
- **Explicit Euler**: simple but limited stability.  
- **Implicit Euler**: robust, stable even for stiff problems.  
- **Implicit Midpoint**: combines stability with improved accuracy.  
