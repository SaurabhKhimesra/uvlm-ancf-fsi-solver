# Fluid–Structure Interaction Analysis using FEM + UVLM  
### **Strong coupling, explicit added‑mass compensation, and energy‑consistent diagnostics**  
**Authors:** Saurabh Khimesra · Ashish Hol  

<p align="center">
  <a href="https://youtu.be/YLkvCXEkd9A">▶️ Demo video</a> ·
  <a href="#quickstart">Quickstart</a> ·
  <a href="#theory-in-3-minutes">Theory</a> ·
  <a href="#reproducing-the-paper-results">Reproduce paper results</a> ·
  <a href="#citation">Citation</a>
</p>

---

## What this repo is 🧠

This repository implements a **time‑domain post‑flutter FSI solver** for thin plates (“flags / flapping sheets”) by coupling:

- **Structure:** geometrically‑nonlinear **ANCF thin‑shell FEM**, time‑integrated with Newmark‑β / generalized‑α
- **Fluid:** **Unsteady Vortex Lattice Method (UVLM)** with wake shedding and convection
- **Coupling:** a **strongly‑coupled staggered partitioned scheme** that stays stable at high added‑mass using  
  **(i)** explicit added‑mass compensation inside the structural solve and **(ii)** Aitken Δ² dynamic relaxation

The accompanying paper (PDF included in this repo) focuses on the numerical stability problem known as **artificial added‑mass instability** and shows how strong coupling + compensation reduces non‑physical energy injection.

**Supported configurations**
- **Single sheet** clamped or pinned at the leading edge  
- **Double sheets** separated by nondimensional spacing \(D^\*\) (interaction effects + wake coupling)

---

## Table of contents

1. [Quickstart](#quickstart)  
2. [Repository layout](#repository-layout)  
3. [Dependencies](#dependencies--toolboxes)  
4. [Running the simulation (GUI workflow)](#running-the-simulation-gui-workflow)  
5. [Key parameters](#key-parameters)  
6. [Theory in 3 minutes](#theory-in-3-minutes)  
7. [Strong coupling algorithm](#strong-coupling-algorithm)  
8. [Coupling diagnostics: energy consistency](#coupling-diagnostics-energy-consistency)  
9. [Reproducing the paper results](#reproducing-the-paper-results)  
10. [Visualization outputs](#visualization-outputs)  
11. [Troubleshooting](#troubleshooting)  
12. [Citation](#citation)  
13. [References used in the paper](#references-used-in-the-paper)  

---

## Quickstart

1. **Clone / download** the repository  
2. Open MATLAB and set the project root as the current directory  
3. Go to either:
   - `single_sheet/` or  
   - `double_sheets/`
4. Add required toolboxes (see [Dependencies](#dependencies--toolboxes))  
5. Open the GUI:
   - open `GUI.fig`
6. Click:
   - **Parameters** → edit
   - **exe** → run solver
   - **plot** → generate figures/movies  
7. Outputs appear under:
   - `./XXXX/save/` (where `XXXX` is `single_sheet` or `double_sheets`)

> If you only want a quick “does it run?” test: keep defaults in `./XXXX/save/param_setting.m` and run ~1–2 seconds of nondimensional time first (smaller `End_Time`) to verify your setup.

---

## Repository layout

```
├─double_sheets
│  ├─cores
│  │  ├─functions
│  │  │  ├─fluid               # UVLM: influence matrix, induced velocity, wake update, force mapping
│  │  │  └─structure           # ANCF shell routines, internal forces, constraints, time stepping helpers
│  │  └─solver
│  │      ├─fluid              # UVLM solve loop (Γ, wake, pressure → nodal forces)
│  │      └─structure           # structural solve loop (Newmark/generalized-α; compensation)
│  ├─save
│  │  └─fig
│  │      └─modes               # saved mode snapshots / modal plots
│  └─ToolBoxes                  # external MATLAB add-ons (see below)
└─single_sheet
   ├─functions
   │  ├─fluid
   │  └─structure
   ├─save
   │  └─fig
   │      └─modes
   ├─solver
   │  ├─fluid
   │  └─structure
   └─ToolBoxes
```

---

## Dependencies / ToolBoxes

The project uses several MATLAB Central submissions for meshing, sparse indexing helpers, and plotting / video export.

### Required (analysis)
- **Meshing a plate using four noded elements** (KSSV)  
  https://www.mathworks.com/matlabcentral/fileexchange/33731-meshing-a-plate-using-four-noded-elements
- **Sparse sub access** (Bruno Luong)  
  https://www.mathworks.com/matlabcentral/fileexchange/23488-sparse-sub-access
- **Vectorized Multi‑Dimensional Matrix Multiplication** (Darin Koblick)  
  https://www.mathworks.com/matlabcentral/fileexchange/47092-vectorized-multi-dimensional-matrix-multiplication

### Recommended (plotting / movies)
- **TriStream** (Matthew Wolinsky)  
  https://www.mathworks.com/matlabcentral/fileexchange/11278-tristream
- **mmwrite** (Micah Richert)  
  https://www.mathworks.com/matlabcentral/fileexchange/15881-mmwrite
- **Quiver 5** (Bertrand Dano)  
  https://www.mathworks.com/matlabcentral/fileexchange/22351-quiver-5

### Add toolbox paths
Edit `add_pathes.m` under `./XXXX/ToolBoxes/` to include the folders you installed, e.g.

```matlab
addpath ./ToolBoxes/<ToolboxName>;
```

### TriStream patch (needed for modern MATLAB versions)
Some MATLAB releases require small edits for stream plotting.

```diff
Line 49:
- x=x(:)'; y=y(:)'; x0=x0(:)'; y0=y0(:)'; u=u(:)'; v=v(:)';
+ x=x(:)'; y=y(:)'; x0=x0(:)'; y0=y0(:)'; u=double(u(:)'); v=double(v(:)');
```

```diff
Line 63:
- TRI = tsearch(x,y,tri',Xbeg,Ybeg);
+ TRI = tsearchn([x.' y.'], tri',[Xbeg.' Ybeg.']);
```

---

## Running the simulation (GUI workflow)

1. Open `GUI.fig` in MATLAB  
2. Click **Parameters** and set:
   - mesh resolution \(N_x \times N_y\)
   - nondimensional parameters \(M^\*\), \(U^\*\), aspect ratio \(H^\*\), thickness \(h^\*\), spacing \(D^\*\) (double case)
   - damping, integrator settings, wake settings
3. Click **exe** to run the time marching + inner coupling loop  
4. Click **plot** to generate the saved figures / movies  
5. Check outputs at `./XXXX/save/` (time histories, snapshots, wake plots, videos)

---

## Key parameters

Defaults / typical analysis settings are stored in:
- `./XXXX/save/param_setting.m`

```matlab
End_Time    = 20;      % Nondimensional analysis time [-]
d_t         = 1.0e-3;  % Nondimensional time step [-]
core_num    = 6;       % CPU core count [-]
speed_check = 0;       % 1: ON, 0: OFF [-]
alpha_v     = 0.5;     % 1: implicit solver, 0: explicit solver [-]

Ma          = 1.0;     % Mass ratio M* [-]
Ua          = 15.0;    % Nondimensional flow velocity U* [-]
theta_a_vec = 0e-1*[0 10]; % Material damping [-]
```

### Nondimensional groups (as defined in the paper)

Let \(L\) be plate length, \(H\) width, \(h\) thickness, \(U_\infty\) free‑stream speed, \(E\) Young’s modulus, and \(I = H h^3/12\).

Mass ratio:
\[
M^\* := \frac{1}{\mu} = \frac{\rho_f L}{\rho_m h}.
\]

Nondimensional flow velocity:
\[
U^\* := \sqrt{\frac{\mu}{\eta}}
      = \sqrt{\frac{\rho_m h L^2 H U_\infty^2}{E I}},
\qquad
I = \frac{H h^3}{12}.
\]

Additional geometry:
- aspect ratio \(H^\* = H/L\)  
- thickness \(h^\* = h/L\)  
- double‑sheet spacing \(D^\* = D/L\)

---

## Theory in 3 minutes

### Governing pieces

**Fluid pressure (inviscid, unsteady Bernoulli):**
\[
p - p_\infty = -\rho_f\left(\frac{\partial\phi}{\partial t} + \frac{1}{2}\|\nabla\phi\|^2\right).
\]

**Structural dynamics (semi‑discrete after FE):**
\[
M\ddot q + C\dot q + f_\text{int}(q) = f_\text{ext}(q,\dot q,t) + f_f(q,\dot q,\Gamma,t).
\]

- \(q\) collects nodal ANCF shell coordinates  
- \(f_f\) is the UVLM aerodynamic load mapped to structural nodes  
- Nonlinearity enters through geometry and internal forces \(f_\text{int}(q)\)

### UVLM in one paragraph

UVLM models inviscid, incompressible flow using bound vortex rings on the plate and a shed wake. At each sub‑iteration (for fixed geometry and wake), the no‑penetration condition at collocation points gives a linear system:

\[
(\mathbf u_\infty + \mathbf u_\text{ind}(\Gamma,\Gamma^w))\cdot \mathbf n = \mathbf v_s\cdot \mathbf n
\quad \Rightarrow \quad
A(q)\Gamma = b(q,\dot q,\Gamma^w).
\]

Wake points convect via:
\[
\mathbf x^w_{n+1} = \mathbf x^w_n + \Delta t\,\mathbf u(\mathbf x^w_n).
\]

---

## Strong coupling algorithm

The coupled problem at each time step can be written abstractly as a fixed‑point residual:

\[
R(q_{n+1}) = 0,
\]
where \(R\) enforces structural dynamics **and** load consistency with the UVLM state.

### Why loose coupling can blow up

In incompressible FSI with high fluid density / low structural mass, the fluid load depends strongly on structural acceleration (added‑mass effect). A naive staggered update can inject a spurious force component and destabilize the scheme (artificial added‑mass instability).

### Core idea: explicit added‑mass compensation inside the structural solve

Approximate the acceleration dependence of the aerodynamic force:
\[
f_f \approx f_f^{(0)} - M_a \ddot q.
\]

Substitute into the structural equation to obtain an “effective inertia” system solved in each sub‑iteration:
\[
(M + M_a)\ddot q + C\dot q + f_\text{int}(q) = f_\text{ext} + f_f^{(0)}.
\]

This shifts the stiff coupling contribution to the implicit side, improving robustness at large \(M^\*\).

### Aitken Δ² dynamic relaxation (fixed‑point accelerator)

Let \(r^{(k)}\) be an interface residual at sub‑iteration \(k\). A relaxed update is:
\[
q^{(k+1)} = q^{(k)} + \omega^{(k)}\Delta q^{(k)}.
\]

Aitken updates \(\omega\) using:
\[
\omega^{(k)} = -\omega^{(k-1)}\,
\frac{(r^{(k)}-r^{(k-1)})^\top r^{(k-1)}}
{\|r^{(k)}-r^{(k-1)}\|^2},
\]
with clipping/safeguards to keep \(\omega\) bounded.

### “Upgrade path” for research‑grade sweeps (from the paper)

For large parameter sweeps (mapping bifurcations / LCO regimes), the paper recommends:
- **Jacobian‑informed interface preconditioning:** update \(M_a\) intermittently (every \(k\) steps or when residual stagnates) rather than every sub‑iteration  
- **Wake continuation / warm‑starting:** extrapolate wake + circulation to reduce inner iterations and improve phase accuracy in LCO regimes

---

## Coupling diagnostics: energy consistency

A practical indicator of coupling fidelity in post‑flutter problems is an energy/work balance check.

Define:
\[
T = \frac{1}{2}\dot q^\top M \dot q,\qquad
V = \Pi(q),\qquad
P_f = f_f^\top \dot q.
\]

Over a time step, the change in mechanical energy should match the work input (minus damping dissipation), up to numerical error.

**Interpretation**
- **Loose coupling** often shows non‑physical energy injection when added‑mass influence is large  
- **Strong coupling + compensation** should reduce those artifacts (see paper Fig. 5 trend)

---

## Reproducing the paper results

The paper demonstrates representative post‑flutter outcomes and visualization (single vs double sheets). A practical workflow:

1. **Start with a stable baseline**
   - moderate \(U^\*\) and \(M^\*\)
   - coarse mesh \(N_x \times N_y\) (fast iteration)
2. **Identify flutter onset region**
   - increase \(U^\*\) gradually
   - monitor tip displacement, dominant frequency, and energy balance
3. **Enter LCO regime**
   - once oscillations persist, refine:
     - time step \(\Delta t\)
     - mesh resolution
     - wake length / truncation
4. **Double sheet studies**
   - vary \(D^\*\) and observe:
     - in‑phase / out‑of‑phase modes
     - wake interference patterns
     - mode switching

### Initial disturbance used to trigger post‑flutter motion

The paper uses a smooth start‑up disturbance of the form:
\[
q_\text{in}(t) = A\sin\left(\frac{\pi t}{t_0}\right)\,\mathbb{I}_{t<t_0},
\]
to avoid numerical shocks and break trivial equilibrium.

---

## Visualization outputs

Typical saved outputs include:
- time histories (tip displacement, selected DOFs)
- wake geometry snapshots / 3D surfaces
- induced‑velocity field visualizations (UVLM output)
- mode snapshots for single and double sheet cases
- movies (if mmwrite installed)

The paper emphasizes qualitative visualization:
- instantaneous streamlines (not streaklines) around sheets  
- wake structure evolution and phase‑lagged loading  
- snapshots of spanwise deformation for double sheets  

---

## Troubleshooting

### 1) Solver diverges / NaNs during coupling
Try, in this order:
- reduce `d_t` (time step)
- increase maximum inner sub‑iterations (coupling iterations)
- enable / strengthen relaxation safeguards (clip \(\omega\) tighter)
- coarsen wake (shorter wake, fewer panels) while debugging
- verify your boundary conditions (clamped vs pinned constraints)

### 2) Wake plots look “explosive”
- reduce time step and/or wake convection velocity options  
- truncate / coarsen far‑wake panels  
- check Kutta condition enforcement and wake roll‑up regularization (if present)

### 3) Streamline / plotting errors (TriStream)
Apply the TriStream patch in the [Dependencies](#tristream-patch-needed-for-modern-matlab-versions) section.

### 4) MATLAB version / compatibility
- the original toolbox set targets broad MATLAB compatibility; if you are on a newer release, the TriStream/tsearch patch is commonly required  
- if `parfor`/parallel features are used, install the Parallel Computing Toolbox (optional)

---

## Citation

**Paper title (PDF in this repo):**  
**Strongly Coupled FEM–UVLM Simulation of Post‑Flutter Limit‑Cycle Oscillations in Thin Flapping Plates with Explicit Added‑Mass Compensation and Energy‑Consistent Coupling**

**BibTeX**
```bibtex
@misc{KhimesraHol_FEMUVLM_StrongCoupling_2025,
  title        = {Strongly Coupled FEM--UVLM Simulation of Post-Flutter Limit-Cycle Oscillations in Thin Flapping Plates with Explicit Added-Mass Compensation and Energy-Consistent Coupling},
  author       = {Khimesra, Saurabh and Hol, Ashish},
  note         = {Manuscript (PDF included in the repository)},
  year         = {2025}
}
```

---

## References used in the paper

1. A. Yamano et al., *Mechanical Engineering Journal* 8(1), 2021. doi:10.1299/mej.20-00459  
2. A. Yamano and M. Chiba, *Int. J. Structural Stability and Dynamics* 22(14), 2022. doi:10.1142/S0219455422501632  
3. A. Yamano et al., *Journal of Sound and Vibration* 478, 2020. doi:10.1016/j.jsv.2020.115359  
4. A. Yamano et al., ICCFD12 paper, 2024 (online PDF link in manuscript).  
5. C. Förster et al., *Comput. Methods Appl. Mech. Eng.* 196(7), 2007. doi:10.1016/j.cma.2006.09.002  
6. P. Causin et al., *Comput. Methods Appl. Mech. Eng.* 194, 2005. doi:10.1016/j.cma.2004.12.005  
7. U. Küttler and W. A. Wall, *Computational Mechanics* 43(1), 2008. doi:10.1007/s00466-008-0255-5  
8. J. Katz and A. Plotkin, *Low‑Speed Aerodynamics*, 2nd ed., Cambridge Univ. Press, 2001.  
9. R. Murua et al., *Progress in Aerospace Sciences* 55, 2012. doi:10.1016/j.paerosci.2012.06.001  
10. A. A. Shabana, *Computational Continuum Mechanics*, Cambridge Univ. Press, 2008.  
11. N. M. Newmark, *J. Eng. Mech. Div., ASCE* 85(EM3), 1959. doi:10.1061/JMCEA3.0000098  
12. J. Chung and G. M. Hulbert, *J. Appl. Mech.* 60(2), 1993. doi:10.1115/1.2900803  
13. M. Chen et al., *Journal of Fluids and Structures* 45, 2014. doi:10.1016/j.jfluidstructs.2013.11.020  

---

## License
See `LICENSE` for licensing terms.  

## Contributing
Issues and pull requests are welcome. Please include:
- MATLAB version + OS  
- the `param_setting.m` you used  
- minimal reproduction steps and (if possible) figures/movies  

