# B-spline Structured Growth Transition Matrix (GTM)

This repository contains the MATLAB implementation developed for **Paper B of my PhD thesis**, where I introduce a B-spline–based Growth Transition Matrix (GTM) model for size-structured fish population dynamics.  
The model provides a flexible and data-driven way to estimate growth transitions and integrate them into a deterministic population projection framework.

The code here reproduces the analyses, model outputs, and figures included in the paper.

---

## 📘 Overview

The methodology includes:

- B-spline estimation of individual growth increments  
- Construction of a probabilistic Growth Transition Matrix (GTM)  
- Models for Normal, Gamma, and Log-normal increment distributions  
- Length-dependent mortality and maturation models  
- A deterministic simulation model to project population structure  
- Tools for comparison of different GTM formulations  

---

## 🚀 How to Run the Code

1. Open MATLAB and go to the project folder:
    ```matlab
    cd /path/to/Bspline-structured-gtm
    ```

2. Initialize the environment:
    ```matlab
    startup
    ```

3. Run the main model (Paper C results):
    ```matlab
    Main
    ```

4. For model comparison (Appendix A) run:
    ```matlab
    Main_comparison
    ```

---

---

## 📝 License
Released under the **MIT License**. You are free to use, modify, and distribute this code with attribution.


