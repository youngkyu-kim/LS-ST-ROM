# Space–Time Reduced Order Model for Linear Dynamical Systems

## 📌 Introduction  
This repository contains the source code for the paper **"Efficient Space–Time Reduced Order Model for Linear Dynamical Systems in Python Using Less than 120 Lines of Code"**, published in *Mathematics* **2021**, 9, 1690 ([DOI: 10.3390/math9141690](https://doi.org/10.3390/math9141690)).  

The code implements a **Space–Time Reduced Order Model (ROM)** to accelerate the simulation of **linear dynamical systems**, leveraging **Proper Orthogonal Decomposition (POD)** with **Galerkin** and **Least-Squares Petrov–Galerkin (LSPG)** projections.

## 🛠 Features  
- Implements both **Galerkin** and **LSPG** projections for space–time ROMs.  
- Handles **diffusion and convection-diffusion** equations in **2D**.  
- Uses **Python** for efficient numerical implementation (less than **120 lines of code** per example).

## 📂 Repository Structure  
```
┌── Scripts/                           # Example scripts for running ROM simulations  
│   ├── Diff_Source_Galerkin.py        # Galerkin ROM for 2D Linear Diffusion Equation with Source Term
│   ├── Diff_Source_LSPG.py            # LSPG ROM for 2D Linear Diffusion Equation with Source Term
│   ├── Conv_Diff_Galerkin.py          # Galerkin ROM for 2D Linear Convection Diffusion Equation
│   ├── Conv_Diff_LSPG.py              # LSPG ROM for 2D Linear Convection Diffusion Equation
│   ├── Conv_Diff_Source_Galerkin.py   # Galerkin ROM for 2D Linear Convection Diffusion Equation with Source Term
│   ├── Conv_Diff_Source_LSPG.py       # LSPG ROM for 2D Linear Convection Diffusion Equation with Source Term
└── README.md                          # This file  
```

## 📖 Citation  
If you use this code in your research, please cite:  
```bibtex
@Article{math9141690,
AUTHOR = {Kim, Youngkyu and Wang, Karen and Choi, Youngsoo},
TITLE = {Efficient Space–Time Reduced Order Model for Linear Dynamical Systems in Python Using Less than 120 Lines of Code},
JOURNAL = {Mathematics},
VOLUME = {9},
YEAR = {2021},
NUMBER = {14},
ARTICLE-NUMBER = {1690},
URL = {https://www.mdpi.com/2227-7390/9/14/1690},
ISSN = {2227-7390},
DOI = {10.3390/math9141690}
}
```

## 📜 License  
This project is licensed under the **Creative Commons Attribution 4.0 (CC BY 4.0)** license.
