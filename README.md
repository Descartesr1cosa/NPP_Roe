# NPP_Roe

A 1D solver for ideal MHD Riemann problems using our new scheme: **Numerical Path Preserving (NPP) Roe scheme**.

---

## Build

This project uses **CMake**.

```bash
mkdir -p build
cd build
cmake ..
make -j
```

Copy the executable `NPP_Roe` to `./0_run`:

```bash
cp ./NPP_Roe ../0_run
```

---

## Run

```bash
cd ../0_run
./NPP_Roe
```

---

## Setup

`NPP_Roe` requires an `input.txt` to specify the simulation setup. Example:

```text
Number_of_Points:              1000
Length_of_Simulation:          3.0
Bx:                            0.7746
CFL:                           0.5
total_time:                    0.3
gamma:                         1.6667
if_NPP(0.or.1)                 0

===========Initial_Left==================
rho--------u-------------v-----------w-----------------By--------------Bz-----------------p
1          0             0.0         0.0              0.7746          0.0             0.5

==========Initial_Right==================
rho--------u------------v-----------w------------------By--------------Bz-----------------p
0.2        0            0.0         0.0               -0.7746         0.1             0.12
```

- `if_NPP` controls the scheme:
  - `0`: Classical Roe scheme
  - `1`: NPP Roe scheme
- Warning: `Bx` should not be too small or zero (not tested by the developer).

---

## Post-processing

- Three test cases are provided at the end of `./0_run/input.txt`.
- After running `NPP_Roe`, the solution will be written to `./output.txt`.
- Copy `output.txt` to the corresponding folder/file for **Example1**, **Example2**, or **Example3**.
- A Tecplot layout file is provided.

**Note**
- Post-processing of **Example1** and **Example2** requires the corresponding **exact MHD Riemann solutions**, which are also provided.
- All provided exact solutions are **REGULAR solutions** (without any intermediate shock).

---

## References

- Exact ideal magnetohydrodynamic Riemann solutions considering the strength of intermediate shocks, DOI: 10.1063/5.0185483  
- Numerical path preserving Godunov schemes for hyperbolic systems, DOI: 10.1016/j.jcp.2023.112297  
- Numerical Path Preserving Roe Scheme for Ideal MHD Riemann Problem: Complete Elimination of Pseudo-Convergence  

---

## Contact

- 1905xuke@buaa.edu.cn
