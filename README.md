# MATLAB Numerical Experiments

This folder contains reproducible MATLAB code for the numerical comparisons in
Section 6 of the manuscript.

Run all experiments from MATLAB with:

```matlab
cd('<paper-folder>/matlab')
run_all_experiments
```

For submission as a single code file, the same computations are also collected
in `paper_numerical_experiments_all_in_one.m`:

```matlab
cd('<paper-folder>/matlab')
results = paper_numerical_experiments_all_in_one;
```

Outputs are written under `matlab/output/`:

- `figures/`: regenerated PNG figures.
- `tables/`: method summaries, pointwise error curves, and chart-condition tables.
- `data/`: MAT files containing the corresponding MATLAB structures.

By default, `run_all_experiments.m` does not overwrite the PNG files referenced
by `main.tex` in the project root.  To overwrite those paper figures, set
`cfg.writePaperFigures = true` in `run_all_experiments.m`.

Implemented comparisons:

- Riemannian normal coordinates with direct Hermite interpolation.
- Raw fixed-block affine chart with direct Hermite interpolation.
- Greedy maximum-volume row-selected affine chart with direct Hermite interpolation.
- CCC chart with direct Hermite interpolation.
- The same chart families with CV+A Hermite interpolation.
- V+A and direct Lagrange utility routines are included for follow-up ablations.
