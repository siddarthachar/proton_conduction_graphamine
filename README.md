# Proton Transport on Graphamine: A Deep-Learning Potential Study

[![DOI](https://img.shields.io/badge/DOI-placeholder-blue)](https://doi.org/10.0000/placeholder)
[![DOI](https://img.shields.io/badge/DOI-placeholder-orange)](https://doi.org/10.0000/placeholder)

## Authors
- Lakshmi Y. Ananthabhotla — Department of Chemical & Petroleum Engineering, University of Pittsburgh, Pittsburgh, PA, 15261, USA
- Siddarth K. Achar — Computational Modeling & Simulation Program, University of Pittsburgh, Pittsburgh, PA, 15260, USA *(Current address: Pritzker School of Molecular Engineering, University of Chicago, Chicago, Illinois 60637, United States)*
- J. Karl Johnson — Department of Chemical & Petroleum Engineering, University of Pittsburgh, Pittsburgh, PA, 15261, USA ([karlj@pitt.edu](mailto:karlj@pitt.edu))

## Abstract
The performance of proton-exchange membrane fuel cells is critically dependent on the conduction of protons. Conventional proton exchange membranes employ materials, such as Nafion, that conduct protons only when properly hydrated. If the relative humidity is too low or too high, the fuel cell will cease to operate. This limitation highlights the need to develop new materials that can rapidly conduct protons without the need to regulate hydration. We present detailed atomistic simulations predicting that graphamine, which is an aminated graphane, conducts protons anhydrously with a very low diffusion barrier compared to existing materials. We have constructed a deep learning framework tailored to modeling graphamine, enabling us to fully characterize and evaluate proton conduction within this material. The trained deep-learning potential is computationally economical and has near-density functional theory accuracy. We used our deep-learning potential to calculate the proton diffusion coefficients at different temperatures and to estimate the activation energy barrier for proton diffusion and found a very low barrier of 63 meV. We estimate the proton conductivity of graphamine to be 1322 mS/cm at 300 K. We show that protons hop along Grotthuss chains containing several amine groups and that the multi-directional hydrogen bonding network intrinsic in graphamine is responsible for the fast conduction of protons.

## Table of Contents
- [Project Overview](#project-overview)
- [Repository Layout](#repository-layout)
- [Environment Setup](#environment-setup)
- [Reproducing Key Workflows](#reproducing-key-workflows)
  - [1. Proton Center-of-Excess Charge Tracking](#1-proton-center-of-excess-charge-tracking)
  - [2. Proton-Conduction Trajectory Unwrapping](#2-proton-conduction-trajectory-unwrapping)
  - [3. Rotational Perturbations for Training Structures](#3-rotational-perturbations-for-training-structures)
- [Data Files](#data-files)
- [Extending the Workflow](#extending-the-workflow)
- [Support](#support)
- [Citation](#citation)

## Project Overview
This repository supports the manuscript *Proton Transport on Graphamine: A Deep-Learning Potential Study*. The code base bundles utilities that were used to construct, validate, and deploy deep learning interatomic potentials for proton migration in graphamine. The workflows draw inspiration from the [reactive_active_learning](https://github.com/siddarthachar/reactive_active_learning) project and adapt its organization to the graphamine-specific simulations reported in the paper.

The key objectives of this repository are to:
- Provide scripts for post-processing molecular dynamics (MD) trajectories of graphamine proton transport.
- Offer tools to generate structural perturbations that increase the diversity of training data for deep-potential models.
- Archive representative structural data and the trained neural network potential that delivers near-DFT accuracy at reduced computational cost.

## Repository Layout
```
README.md                     Project documentation.
cec.py                        Center-of-excess-charge analysis for proton tracking.
gamine_DP.pb                  Trained deep-learning potential in DeepMD format.
gamine_np.poscar              Reference neutral protonated POSCAR structure.
gamine_p.poscar               Reference protonated POSCAR structure.
rotate.py                     Generates rotated NH₂ configurations for sampling.
xdatcar2unwrappedxyz.py       Converts XDATCAR trajectories to unwrapped XYZ.
```

## Environment Setup
All Python utilities were written for Python 3.8+ and leverage the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/). A minimal environment can be created with `conda` or `pip`:

```bash
conda create -n graphamine python=3.10
conda activate graphamine
pip install ase numpy matplotlib
```

Some scripts expect additional tooling (e.g., LAMMPS for generating trajectories) and external helper scripts referenced via absolute paths in the original workflow. When porting the code to a new environment, adapt those paths to your local setup.

## Reproducing Key Workflows

### 1. Proton Center-of-Excess Charge Tracking
`cec.py` computes the nitrogen-centered excess charge (CEC) trajectory of the proton as it hops between amine sites.

**Inputs**
- `dumptraj`: LAMMPS dump or VASP XDATCAR trajectory containing scaled proton positions.
- `dumptraj_tag`: Either `"LAMMPS"` or `"XDATCAR"` to control parsing.
- `noAtoms`: Number of atoms in the simulation cell.
- `xdatcar_prepend`: Header string for the generated `XDATCAR_CEC` file.

**Usage sketch**
```python
from cec import CEC
cec = CEC(
    working_dir='path/to/workdir',
    dumptraj='trajectory.lammpstrj',
    dumptraj_tag='LAMMPS',
    noAtoms=number_of_atoms,
    xdatcar_prepend=open('XDATCAR_header').read(),
    timestep=0.25,
    rcut=1.3,
)
cec.CECmain('cec_indices.txt')
```
The script identifies the nitrogen hosting the proton at each step, writes a proton-hopping `XDATCAR_CEC`, and stores the time series of nitrogen indices for downstream analysis.

### 2. Proton-Conduction Trajectory Unwrapping
`xdatcar2unwrappedxyz.py` converts the `XDATCAR_CEC` output into an unwrapped `coordunwrapped_CEC.xyz` trajectory suitable for mean-squared displacement (MSD) calculations. Run the script with:

```bash
python xdatcar2unwrappedxyz.py -i XDATCAR_CEC
```

The current implementation references a `POSCAR_prod` file through an absolute path. Update the path to point to your local POSCAR if needed. The generated XYZ trajectory can then be processed with external MSD utilities (see `cec.py` for the expected command).

### 3. Rotational Perturbations for Training Structures
`rotate.py` generates ensembles of rotated NH₂ configurations starting from a pristine `POSCAR`. These perturbations enrich the structural diversity of the training set for the deep potential.

**Key steps performed by the script**
1. Identify hydrogen and nitrogen indices from the input `POSCAR`.
2. Construct neighbor lists to pair each nitrogen with its two nearest hydrogens.
3. Create `rotate_no` copies of the base structure and rotate the hydrogen pairs by randomized angles around the nitrogen centers.

To execute:
```bash
python rotate.py
```
The script writes `POSCAR-0`, `POSCAR-1`, … files containing the rotated structures. Adjust `rotate_no` at the top of the script to control how many variants are generated and revise hard-coded directories before running in a new environment.

## Data Files
- **`gamine_DP.pb`** — Trained Deep Potential (DeepMD) model capturing graphamine proton dynamics. Load this with DeepMD-kit or LAMMPS (`pair_style deepmd`).
- **`gamine_np.poscar`** — Reference neutral protonated graphamine structure used in training and validation.
- **`gamine_p.poscar`** — Protonated configuration highlighting the migration pathway used to seed MD simulations.

## Extending the Workflow
1. **Swap in new trajectories:** Point `cec.py` to alternative LAMMPS or VASP outputs to analyze additional proton transport simulations.
2. **Customize neighbor heuristics:** Modify the `rcut` parameter in `CEC` to tune proton-host detection for systems with different bonding distances.
3. **Integrate with Deep Potential MD:** Use `gamine_DP.pb` with DeepMD-enabled MD engines to run large-scale proton conduction simulations at reduced cost.
4. **Automate MSD calculations:** The `calc_msd` stub in `cec.py` illustrates how to bridge into a mean-squared-displacement utility. Replace the placeholder paths with your preferred analysis scripts.

## Support
For questions about the workflows or to request additional data, please contact [karlj@pitt.edu](mailto:karlj@pitt.edu).

## Citation
If you use this repository in your research, please cite:

> Lakshmi Y. Ananthabhotla, Siddarth K. Achar, and J. Karl Johnson. *Proton Transport on Graphamine: A Deep-Learning Potential Study.* University of Pittsburgh.

Placeholders for DOIs are provided above; update them when an official DOI becomes available.
