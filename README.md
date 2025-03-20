# Efficient Space–Time Reduced Order Model for Linear Dynamical Systems

![License](https://img.shields.io/badge/license-CC%20BY%204.0-blue.svg)  
[![DOI](https://zenodo.org/badge/DOI/10.3390/math9141690.svg)](https://doi.org/10.3390/math9141690)

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
├── README.md                          # This file  
└── LICENSE                            # CC-BY-4.0 license  
```

## 🚀 Quick Start  
### Prerequisites  
Ensure you have **Python 3.7+** installed with the following dependencies:  
```bash
pip install numpy scipy matplotlib
```

### Running an Example  
To run the **2D diffusion equation** example:  
```bash
python examples/diffusion.py
```
For the **2D convection-diffusion** equation:  
```bash
python examples/convection_diffusion.py
```

## 📊 Results & Performance  
The ROM implementation achieves speed-ups of **O(10-1000x** compared to full-order models (FOMs) while maintaining relative errors below **0.01%** in most cases.  

Sample results from the paper:  
![Performance Comparison](results/performance_plot.png)

## 📖 Citation  
If you use this code in your research, please cite:  
```bibtex
@article{kim2021efficient,
  author = {Kim, Youngkyu and Wang, Karen and Choi, Youngsoo},
  title = {Efficient Space–Time Reduced Order Model for Linear Dynamical Systems in Python Using Less than 120 Lines of Code},
  journal = {Mathematics},
  volume = {9},
  number = {14},
  pages = {1690},
  year = {2021},
  publisher = {MDPI},
  doi = {10.3390/math9141690}
}
```

## 📜 License  
This project is licensed under the **Creative Commons Attribution 4.0 (CC BY 4.0)** license. See the [LICENSE](LICENSE) file for details.
