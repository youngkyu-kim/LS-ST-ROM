# Efficient Space–Time Reduced Order Model for Linear Dynamical Systems

![License](https://img.shields.io/badge/license-CC%20BY%204.0-blue.svg)  
[![DOI](https://zenodo.org/badge/DOI/10.3390/math9141690.svg)](https://doi.org/10.3390/math9141690)

## 📌 Introduction  
This repository contains the source code for the paper **"Efficient Space–Time Reduced Order Model for Linear Dynamical Systems in Python Using Less than 120 Lines of Code"**, published in *Mathematics 2021, 9, 1690* ([DOI: 10.3390/math9141690](https://doi.org/10.3390/math9141690)).  

The code implements a **Space–Time Reduced Order Model (ROM)** to accelerate the simulation of **linear dynamical systems**, leveraging **Proper Orthogonal Decomposition (POD)** and **Least-Squares Petrov–Galerkin (LSPG) projection**.

## 🛠 Features  
- Implements both **Galerkin** and **LSPG** projections for space–time ROMs.  
- Handles **diffusion and convection-diffusion** equations in **2D**.  
- Uses **Python** for efficient numerical implementation (less than **120 lines of code** per example).  
- Achieves **significant speed-ups** in solving large-scale linear dynamical systems.  
- Includes **error analysis** and **comparisons** between different ROM formulations.

## 📂 Repository Structure  
```
├── examples/                  # Example scripts for running ROM simulations  
│   ├── diffusion.py           # 2D Diffusion problem implementation  
│   ├── convection_diffusion.py # 2D Convection-diffusion problem implementation  
│   ├── utils.py               # Helper functions for numerical operations  
├── data/                      # Sample datasets used in the paper  
├── results/                   # Outputs and plots from the ROM simulations  
├── README.md                  # This file  
└── LICENSE                    # CC-BY-4.0 license  
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
