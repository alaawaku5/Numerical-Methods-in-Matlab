# Project: Euler and Midpoint Methods for ODEs (Exercise 4)

## 📌 Overview
This project implements and compares three numerical schemes in MATLAB for solving linear initial value problems (IVPs):

- Explicit Euler
- Implicit Euler
- Implicit Midpoint

We apply these methods to two test problems:

### Part (a)
y'(t) = -y(t),   y(0) = 1  
Exact solution: y(t) = exp(-t)

### Part (b)
y'(t) = f1(t) + f2(t) * y(t),   y(0) = 5  
with constants f1(t) = 2.0, f2(t) = -1.0

---

## 🧮 Methods
- **Explicit Euler**  
  y[n+1] = y[n] + h * ( f1(t[n]) + f2(t[n]) * y[n] )

- **Implicit Euler**  
  y[n+1] = ( y[n] + h * f1(t[n+1]) ) / ( 1 - h * f2(t[n+1]) )

- **Implicit Midpoint**  
  y[n+1] = ( y[n] + h * f1(t_m) + (h * f2(t_m)/2) * y[n] ) / ( 1 - (h * f2(t_m)/2) )  
  where t_m = ( t[n] + t[n+1] ) / 2

---

## 📊 Results
- **Part (a):**  
  - Explicit Euler: conditionally stable, less accurate  
  - Implicit Euler: stable but slightly diffusive  
  - Implicit Midpoint: stable and more accurate  

- **Part (b):**  
  - Similar behavior observed with forcing term included  
  - Implicit schemes outperform explicit Euler in terms of stability  

### Example Output (Part a)
![Euler Methods Comparison](results_part_a.png)

### Example Output (Part b)
![Euler Methods Comparison](results_part_b.png)

---

## 📂 Files
- `exercise4.m` → MATLAB code (both parts a and b, with exact solution for part a)  
- `results_part_a.png` → plot comparing exact and numerical solutions for part a  
- `results_part_b.png` → plot comparing explicit, implicit, and midpoint for part b  

---

## ✅ Takeaway
This project highlights how explicit and implicit methods differ in accuracy and stability for linear ODEs:

- Explicit Euler → simple but limited stability  
- Implicit Euler → robust, stable even for stiff problems  
- Implicit Midpoint → combines stability with improved accuracy
