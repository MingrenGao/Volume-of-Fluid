# Volume-of-Fluid

A standalone C++ test program to study **curvature accuracy and grid convergence**
for 2D interfaces represented by **VOF (Volume of Fluid)** on a Cartesian grid.

---

## Analytic Interface Shapes

The code supports two analytic interface shapes used for curvature validation.

### Circle

Constant exact curvature:

$$
\kappa = \frac{1}{R}
$$

where \( R \) is the circle radius.

---

### Cosine

Analytic interface defined as:

$$
y(x) = A - \cos\left(\frac{2\pi n x}{L}\right)
$$

with known analytic curvature.

---

## Curvature Methods Implemented

- **CV**  
  Convolved VOF method

- **RDF**  
  Reconstructed distance function method

- **PARABOLOID**  
  Parabolic reconstruction

- **HF**  
  Height-function method

- **VFM**  
  Volumetric fitting method

- **PC**  
  PLIC-centroidal fitting method

- **PV**  
  PLIC-volumetric fitting method

---

## Convergence Study

The program performs a resolution sweep specified by `amr.n_list` and computes:

- **L2 error norm**
- **L∞ error norm**

Errors are evaluated on interface cells only, using the analytic curvature as reference.

At the **finest resolution**, the following diagnostic files are written:

- `kappa_field.dat` — curvature field dump  
- `phi_field.dat` — optional (for RDF methods)  
- `interface_segments.dat` — reconstructed PLIC segments  

---

## Quickstart

### Build (g++)

```bash
g++ -O3 -std=c++17 -Wall -Wextra -pedantic VOF.cpp -o curv
./curv cur.inp
