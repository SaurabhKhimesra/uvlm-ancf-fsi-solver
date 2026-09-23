# Fluid-Structure Interaction of Thin Flapping Plates (ANCF shell + UVLM)

Time-domain solver for post-flutter limit-cycle oscillations of thin plates ("flags"),
coupling a geometrically-nonlinear ANCF shell FEM to an unsteady vortex lattice method.

**Authors:** Saurabh Khimesra, Ashish Hol
**Demo video:** https://youtu.be/YLkvCXEkd9A

MATLAB, no compiled dependencies. Two cases are included: a single sheet, and two
sheets separated by a nondimensional gap.

---

## Contents

- [What it solves](#what-it-solves)
- [Requirements](#requirements)
- [Repository layout](#repository-layout)
- [Running a case](#running-a-case)
- [Parameters](#parameters)
- [Formulation](#formulation)
- [Outputs](#outputs)
- [Troubleshooting](#troubleshooting)
- [References](#references)
- [License](#license)

---

## What it solves

- **Structure** - thin shell discretised with the absolute nodal coordinate
  formulation (ANCF), 9 DOF per node, 36 per four-node element. Large rotations and
  large displacements are handled without small-angle assumptions.
- **Fluid** - unsteady vortex lattice method. Bound vortex rings on the plate, a shed
  wake that convects with the local induced velocity, and an unsteady Bernoulli
  pressure recovered on each panel.
- **Coupling** - the added-mass part of the fluid load is computed explicitly and
  moved onto the inertia side of the structural equation, which is what keeps the
  scheme stable at mass ratios where a plain staggered exchange diverges. Set with
  `coupling_flag` (single-sheet case).

Time integration is an **Euler predictor-corrector** scheme
(`cores/solver/structure/solve_structure.m`), with `alpha_v` blending the implicit and
explicit updates. `solve_structure_simple_implicit.m` and
`solve_structure_PredictorCorrector_method.m` are kept alongside it as alternatives.

---

## Requirements

MATLAB R2016b or later. Parallel Computing Toolbox is optional -
`maxNumCompThreads(core_num)` is used, not `parfor`.

The solver depends on six MATLAB Central submissions. They are **not** vendored here;
download them and drop each one into the `ToolBoxes` folder under the folder name that
`add_pathes.m` expects:

| Folder name expected | Submission | Used by |
|---|---|---|
| `Meshing_a_plate_using_four_noded_elements` / `Plate_Mesh` | [Meshing a plate using four noded elements](https://www.mathworks.com/matlabcentral/fileexchange/33731) (KSSV) | mesh generation |
| `Sparse_sub_access` / `FEM_sparse` | [Sparse sub access](https://www.mathworks.com/matlabcentral/fileexchange/23488) (B. Luong) | global matrix assembly |
| `Vectorized_Multi-Dimensional_Matrix_Multiplication` / `mntimes` | [Vectorized multi-dimensional matrix multiplication](https://www.mathworks.com/matlabcentral/fileexchange/47092) (D. Koblick) | element-wise kernels |
| `TriStream` | [TriStream](https://www.mathworks.com/matlabcentral/fileexchange/11278) (M. Wolinsky) | streamline plots |
| `quiver5` | [Quiver 5](https://www.mathworks.com/matlabcentral/fileexchange/22351) (B. Dano) | velocity field plots |
| `mmwrite` | [mmwrite](https://www.mathworks.com/matlabcentral/fileexchange/15881) (M. Richert) | movie export |

The two cases were written at different times and expect slightly different folder
names - check `single_sheet/add_pathes.m` and `double_sheets/add_pathes.m` and rename
to match. `single_sheet` additionally lists `mpg_write/src` (alternative movie export)
and a commented-out `lightspeed` entry; neither is required. The last three rows above
are only needed for plotting and movies - the solver itself runs without them.

`ToolBoxes/toolboxes.pdf` lists the same set with screenshots.

---

## Repository layout

```
single_sheet/
  GUI.m, GUI.fig           entry point
  add_pathes.m             addpath list (edit to match your ToolBoxes folder names)
  cores/
    exe.m                  main time-marching script
    initializing.m         default figure/font settings
    plot_data.m            post-processing and figure export
    version_check.m        guards param_setting against the solver version
    reduce_mat_file.m      shrinks saved .mat results
    solver/
      fluid/               panels, wake, induced velocity, fluid force (8 files)
      structure/           shape functions, elements, stiffness, time stepping (9 files)
      solve_energy.m       energy and work-rate balance
      solve_mode.m         modal analysis
      initial_values.m
    functions/
      fluid/               vortex-ring kernels, pressure interpolation
      structure/           strain/curvature derivatives, system residual
    ToolBoxes/             <- third-party submissions go here
  save/
    param_setting.m        all run parameters
    fig/                   exported figures (.fig/.pdf), fig/modes/ for mode shapes

double_sheets/             same structure; ToolBoxes/ sits at the case root, not under cores/
```

---

## Running a case

1. Open MATLAB and `cd` into `single_sheet` or `double_sheets`.
2. Install the toolboxes above and reconcile the names in `add_pathes.m`.
3. Edit `save/param_setting.m`.
4. Run `GUI` from the command window (or open `GUI.fig`), then use **Parameters** ->
   **exe** -> **plot**.

To run headless, skip the GUI and call `exe` directly - it calls `add_pathes` and
`param_setting` itself.

`version_check.m` compares `exe_ver` against `param_ver` and warns on a mismatch, so
keep the two in step if you copy a parameter file between the cases (single sheet is
at version 12.0, double sheets at 1.0 - they are not interchangeable).

For a first run, cut `End_Time` to 1-2 and keep the default mesh. A full
`End_Time = 30` run at `Nx = 15, Ny = 10` takes hours.

---

## Parameters

Everything lives in `save/param_setting.m`. Defaults as committed:

```matlab
%% single_sheet
End_Time      = 30;       % nondimensional analysis time
d_t           = 1.5e-3;   % nondimensional time step
core_num      = 8;        % computational threads
alpha_v       = 0.5;      % 1: implicit, 0: explicit
coupling_flag = 1;        % 1: strong coupling (explicit added mass), 0: weak (staggered)

Ma            = 1.0;      % mass ratio M*
Ua            = 25;       % nondimensional flow velocity U*
Nx            = 15;       % elements along the chord
Ny            = 10;       % elements along the span
thick         = 1e-3;     % nondimensional thickness
Width         = 1.0;      % aspect ratio H* = H/L
mode_num      = 5;
```

```matlab
%% double_sheets
End_Time      = 20;
d_t           = 1.0e-3;
core_num      = 6;
Ua            = 15.0;
theta_a_vec   = 0e-1*[0 10];   % material damping, per sheet
```

The initial disturbance is a half-sine body force applied over the first 0.2 of
nondimensional time,

```matlab
q_in_norm = @(time)( 0.5*sin(pi*time/0.2).*(time < 0.2) );
```

which excites the plate without the numerical shock that an impulsive initial velocity
(`dt_rz_end`) produces.

### Nondimensional groups

With $L$ the plate length, $H$ the width, $h$ the thickness, $U_\infty$ the free-stream
speed, $E$ Young's modulus and $I = H h^3 / 12$:

$$M^* = \frac{\rho_f L}{\rho_m h}, \qquad U^* = \sqrt{\frac{\rho_m h L^2 H U_\infty^2}{E I}}$$

plus the geometric ratios

$$H^* = H/L, \qquad h^* = h/L, \qquad D^* = D/L$$

where $D$ is the gap between the two sheets in the double-sheet case.

---

## Formulation

**Structure.** After discretisation the semi-discrete system is

$$M \ddot{q} + C \dot{q} + f_{int}(q) = f_{ext} + f_f(q, \dot{q}, \Gamma, t)$$

where $q$ collects the ANCF nodal coordinates (position and slope vectors) and
$f_{int}$ carries the geometric nonlinearity through the strain and curvature terms.

**Fluid.** Enforcing no-penetration at the panel collocation points gives, for a frozen
geometry and wake,

$$A(q)\Gamma = b(q, \dot{q}, \Gamma^w)$$

and the wake nodes convect as $x^w_{n+1} = x^w_n + \Delta t \cdot u(x^w_n)$. Pressure comes
from the unsteady Bernoulli equation and is integrated to nodal forces.

**Added mass.** The fluid load depends on the structural acceleration. Writing that part
out explicitly,

$$f_f \approx f_f^{(0)} - M_a \ddot{q}$$

and folding $M_a$ into the left-hand side gives an effective inertia
$(M + M_a)$. In the code this is the `Qf_p_mat_global` block assembled in
`calc_fluid_force.m` and applied in `solve_structure.m`; `calc_fluid_force_strong.m`
and `calc_fluid_force_weak.m` hold the two variants. Without it, runs at $M^* \sim 1$
diverge within a few hundred steps.

**Energy check.** `solve_energy.m` accumulates kinetic energy
$T = \tfrac{1}{2}\dot{q}^\top M \dot{q}$, strain energy, damping dissipation and the
fluid work rate $P_f = f_f^\top \dot{q}$ each step. Plotting the residual of that
balance (`save/fig/work_rate.fig`) is the quickest way to tell whether a run is
physically meaningful or being driven by coupling error.

---

## Outputs

Written to `save/` and `save/fig/`:

- `displacement.fig`, `displacement_mid_span.fig` - tip and mid-span time histories
- `disp_vel_mid_span_phase_plane.fig` - phase portrait, shows the limit cycle
- `snapshot.fig`, `snapshot_mid_span.fig` - deformed shapes through one cycle
- `Velocity_field.fig`, `u/v/Unorm_distribution_0.fig` - UVLM induced velocity field
- `work_rate.fig` - energy balance
- `alpha_vs_CL.fig` - lift coefficient against angle of attack
- `fig/modes/mode_*.fig` - mode shapes from `solve_mode.m`
- `data.wmv` - animation, if `mmwrite` is installed

`reduce_mat_file.m` strips the per-step history out of a saved result when the `.mat`
gets unwieldy.

---

## Troubleshooting

**Diverging or NaN after a few steps.** Check `coupling_flag = 1` first - weak coupling
is not stable at the default mass ratio. Then reduce `d_t`, then coarsen the wake.

**`Undefined function` on the first run.** `add_pathes.m` is adding folders that do not
exist. The `ToolBoxes` folders in this repo hold documentation only; the submissions
have to be downloaded separately (see [Requirements](#requirements)).

**Wake blows up visually.** Usually `d_t` too large relative to the panel size, or a
wake that has been left to grow unbounded over a long run.

**`tsearch` undefined in TriStream.** `tsearch` was removed from MATLAB. Patch
`TriStream.m`:

```diff
 line 49:
-x=x(:)'; y=y(:)'; x0=x0(:)'; y0=y0(:)'; u=u(:)'; v=v(:)';
+x=x(:)'; y=y(:)'; x0=x0(:)'; y0=y0(:)'; u=double(u(:)'); v=double(v(:)');

 line 63:
-TRI = tsearch(x,y,tri',Xbeg,Ybeg);
+TRI = tsearchn([x.' y.'], tri',[Xbeg.' Ybeg.']);
```

This is also what `NOTE.pdf` contains.

**Source encoding.** Most `.m` files carry Shift-JIS comments from the original
working version. They run fine, but set MATLAB's encoding to Shift-JIS, or expect
mojibake in the comments. `.gitattributes` marks `.m` as `-text` so git leaves the
bytes untouched on checkout.

---

## References

1. A. Yamano et al., *Mechanical Engineering Journal* **8**(1), 2021. doi:10.1299/mej.20-00459
2. A. Yamano, M. Chiba, *Int. J. Structural Stability and Dynamics* **22**(14), 2022. doi:10.1142/S0219455422501632
3. A. Yamano et al., *Journal of Sound and Vibration* **478**, 2020. doi:10.1016/j.jsv.2020.115359
4. C. Forster, W. A. Wall, E. Ramm, *Comput. Methods Appl. Mech. Eng.* **196**(7), 2007. doi:10.1016/j.cma.2006.09.002
5. P. Causin, J. F. Gerbeau, F. Nobile, *Comput. Methods Appl. Mech. Eng.* **194**, 2005. doi:10.1016/j.cma.2004.12.005
6. J. Katz, A. Plotkin, *Low-Speed Aerodynamics*, 2nd ed., Cambridge University Press, 2001.
7. R. Murua, R. Palacios, J. M. R. Graham, *Progress in Aerospace Sciences* **55**, 2012. doi:10.1016/j.paerosci.2012.06.001
8. A. A. Shabana, *Computational Continuum Mechanics*, Cambridge University Press, 2008.
9. M. Chen et al., *Journal of Fluids and Structures* **45**, 2014. doi:10.1016/j.jfluidstructs.2013.11.020

---

## License

MIT, see [LICENSE](LICENSE). The MATLAB Central submissions listed under
[Requirements](#requirements) are covered by their own licenses and are not
redistributed here.
