# Numerical Methods for ODEs

Implementation and comparison of several multistep numerical methods for approximating the solution of an ordinary differential equation (ODE), written in Python.

## Problem

The test problem is the initial value problem

$$y'(t) = -y(t) + e^{-t}\cos(t), \qquad y(0) = 0, \qquad t \in [0, 5]$$

whose exact solution is known:

$$y(t) = e^{-t}\sin(t)$$

Knowing the exact solution allows us to measure the error of each method precisely.

## Implemented methods

| Method | Function | Starting procedure |
|---|---|---|
| Adams–Bashforth, 2 steps (AB2) | `adams_bashforth2` | Euler's method |
| Adams–Bashforth, 3 steps (AB3) | `adams_bashforth3` | Midpoint method |
| Adams–Moulton, 3 steps (AM3), implicit | `adams_moulton3` | Runge–Kutta 4 (RK4) |

The implicit equation in AM3 is solved at each step with a fixed-point iteration (tolerance `1e-12`, maximum 200 iterations).

## Experiments

1. **Fixed mesh (N = 40):** prints the maximum error of each method and plots all approximations together with the exact solution.
2. **Several meshes (N = 10, 20, 40, 80, 160, 320):** for each method, prints the maximum error for every mesh and plots the approximations.
3. **Orders of convergence:** estimates the empirical order of each method as `log(e_{N/2} / e_N) / log(2)` when the mesh is refined.

## Requirements

- Python 3
- NumPy
- Matplotlib

```bash
pip install -r requirements.txt
```

## Usage

```bash
python Numeric_methodes_for_ODEs.py
```

The script runs all three experiments, printing the results to the console and opening the corresponding figures. The console output and plot labels are in Spanish.

## Author

Juan Antonio Gaspar Pascual
