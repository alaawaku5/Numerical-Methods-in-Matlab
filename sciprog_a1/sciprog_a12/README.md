
# Project: One-Sided Finite Difference Approximation (sciprog_a12.m)

## 📌 Overview
Implements the one-sided finite difference method in MATLAB to approximate derivatives at x = 0.  
Analyzes error behavior and convergence rates for three test functions:
- exp(x), f'(0) = 1
- sin(x), f'(0) = 1
- x^(3/2), f'(0) = 0

## 🧮 Method


- For smooth functions: first-order accurate.  
- For x^(3/2): reduced convergence due to non-smoothness.  

## 📊 Results
- Verified first-order convergence for exp(x), sin(x).  
- Observed reduced convergence for x^(3/2).  

![Error vs Step Size](results.png)

## 📂 Files
- `sciprog_a12.m` → MATLAB code  
- `results.png` → error plot  

## ✅ Takeaway
Showcases numerical error analysis, convergence study, and MATLAB plotting.
