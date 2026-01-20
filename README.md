# Volume-of-Fluid

A standalone C++ test program to study **curvature accuracy and grid convergence** for 2D interfaces represented by **VOF (Volume of Fluid)** on a Cartesian grid.

It supports two analytic interface shapes:

- **circle**: constant exact curvature  

  $$
  \kappa = \frac{1}{R}
  $$

- **cosine**:  

  $$
  y(x) = A - \cos\left( \frac{2\pi n x}{L} \right)
  $$


Curvature methods implemented:

- **CV**: convolution-smoothed VOF + divergence of normals  
- **RDF_raw / RDF_smooth**: reconstructed distance function (optionally smoothed) + level-set curvature
- **PARABOLOID**: local height-function / parabola reconstruction (paper-aligned)
- **HF**: height-function curvature (cosine uses full-column HF for stability)
- **VFM**: volume-fraction matching (LM fit in local tangent-normal frame)
- **PC**: PLIC-segment centroid parabola fit
- **PV**: PLIC-segment moment-based parabola fit (PV variant)

The program sweeps resolutions (`amr.n_list`), computes **L2** and **L∞** error norms on interface cells, and writes diagnostic dumps at the finest grid:
`kappa_field.dat`, optional `phi_field.dat`, `interface_segments.dat`, and optional `alpha.ppm` / `curvature.ppm`.

---

## Quickstart

### Build (g++)

```bash
g++ -O3 -std=c++17 curvature_convergence_paraboloid_fixed.cpp -o curv_test
