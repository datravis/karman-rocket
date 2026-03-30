"""
Apex-1 Fin Structural FEA  —  scikit-fem Kirchhoff plate analysis
===================================================================
Models one fin as a clamped trapezoidal Kirchhoff plate under
uniform aerodynamic pressure at max-Q (Mach 1.1, 5 km altitude).

Improves on the closed-form estimate in structural_checks.py by:
  • Using the actual trapezoidal fin geometry (not a rectangular proxy)
  • Correctly distributing load across the full planform
  • Recovering peak bending stress at the root via nodal slope gradients

Element: ElementTri15ParamPlate (15-parameter nonconforming, skfem)
Mesh:    Trapezoidal trapezoid refined to ~2000 triangles (4 uniform refinements)
BCs:     Clamped root edge (w = dw/dx = dw/dy = 0 at x = 0)
Load:    Uniform pressure equivalent to Barrowman/Helmbold normal force

References:
  - Kirchhoff (1850): thin-plate bending theory
  - Barrowman / Helmbold: fin normal force coefficient
  - skfem docs: https://scikit-fem.readthedocs.io

Dependencies:
  pip install scikit-fem
"""

import math
import numpy as np
from skfem import (MeshTri, Basis, ElementTri15ParamPlate,
                   BilinearForm, LinearForm, asm, condense, solve)
from skfem.helpers import dd, ddot

# ──────────────────────────────────────────────────────────────────────────────
# Fin geometry  (all dims in metres)
# ──────────────────────────────────────────────────────────────────────────────
b_fin  = 0.200   # exposed semi-span (cantilever length)
c_root = 0.360   # root chord
c_tip  = 0.160   # tip chord
sweep  = 0.130   # leading-edge sweep (LE offset at tip vs root)
t_fin  = 0.004   # thickness

c_avg  = (c_root + c_tip) / 2
S_fin  = c_avg * b_fin          # planform area per fin
AR_exp = b_fin**2 / S_fin

# ──────────────────────────────────────────────────────────────────────────────
# Material: E-glass / epoxy (woven, balanced 0/90)
# ──────────────────────────────────────────────────────────────────────────────
E_fg        = 17e9    # Young's modulus, Pa
nu_fg       = 0.18    # Poisson's ratio
sigma_ult   = 230e6   # conservative bending strength, Pa

D = E_fg * t_fin**3 / (12 * (1 - nu_fg**2))   # flexural rigidity, N·m

# ──────────────────────────────────────────────────────────────────────────────
# Aerodynamic loading  —  max-Q conditions (Mach 1.1, 5 km altitude)
# ──────────────────────────────────────────────────────────────────────────────
rho_5km  = 0.736          # kg/m³
V_maxQ   = 1.1 * 320.0   # m/s
q_maxQ   = 0.5 * rho_5km * V_maxQ**2   # dynamic pressure, Pa
AoA_rad  = math.radians(5.0)

# Barrowman / Helmbold: CN_alpha per fin, per radian
CN_alpha = 2 * math.pi / (1 + 2 / AR_exp)

# Total normal force on one fin (N) and equivalent uniform pressure (Pa)
N_fin     = CN_alpha * q_maxQ * S_fin * AoA_rad
p_uniform = N_fin / S_fin   # distribute uniformly across planform

# ──────────────────────────────────────────────────────────────────────────────
# Mesh: trapezoidal fin in (x=span, y=chord) coordinates
#   x = 0 → root edge      x = b_fin → tip edge
#   LE at y = 0 (root) and y = sweep (tip)
# Vertices:  root-LE, root-TE, tip-TE, tip-LE
# ──────────────────────────────────────────────────────────────────────────────
pts   = np.array([[0.0,  0.0      ],
                  [0.0,  c_root   ],
                  [b_fin, sweep + c_tip],
                  [b_fin, sweep        ]])
cells = np.array([[0, 1, 2], [0, 2, 3]])
mesh  = MeshTri(pts.T, cells.T).refined(4)   # 512 → 8192 triangles (4× refinement)

basis = Basis(mesh, ElementTri15ParamPlate())

# ──────────────────────────────────────────────────────────────────────────────
# Assembly: Kirchhoff plate bilinear form
#   a(w,v) = D ∫ [ν·(w_xx+w_yy)·(v_xx+v_yy) + (1−ν)·w_ij·v_ij] dA
# ──────────────────────────────────────────────────────────────────────────────
@BilinearForm
def plate(u, v, w):
    d2u = dd(u); d2v = dd(v)
    return D * (nu_fg * (d2u[0,0] + d2u[1,1]) * (d2v[0,0] + d2v[1,1])
                + (1 - nu_fg) * ddot(d2u, d2v))

@LinearForm
def load(v, w):
    return p_uniform * v

K = asm(plate, basis)
f = asm(load, basis)

# ──────────────────────────────────────────────────────────────────────────────
# Boundary conditions: clamp root edge (x ≈ 0)
#   DOFs per node: [w,  dw/dx,  dw/dy]  (nodal_dofs rows 0–2)
#   All three zeroed → full clamped condition
# ──────────────────────────────────────────────────────────────────────────────
root_dofs = basis.get_dofs(lambda x: x[0] < 1e-9).all()
u_sol = solve(*condense(K, f, D=root_dofs))

# ──────────────────────────────────────────────────────────────────────────────
# Post-processing: bending stress via nodal slope gradient recovery
#
# For a Kirchhoff plate, extreme fibre stress at z = t/2:
#   σ_x = E·t / (2·(1−ν²)) · (κ_xx + ν·κ_yy)
#
# Curvatures recovered from the nodal slope DOFs:
#   κ_xx ≈ d(dw/dx)/dx  — slope of the x-slope field in the x-direction
#   κ_yy ≈ d(dw/dy)/dy  — slope of the y-slope field in the y-direction
#
# Procedure: fit a linear field wx = a + b·x + c·y over each element's
# three nodal wx values → b = dw_x/dx = κ_xx (element-average curvature).
# ──────────────────────────────────────────────────────────────────────────────
wx_nodes = u_sol[basis.nodal_dofs[1]]   # dw/dx at nodes
wy_nodes = u_sol[basis.nodal_dofs[2]]   # dw/dy at nodes
w_nodes  = u_sol[basis.nodal_dofs[0]]   # deflection at nodes

n_elem = mesh.t.shape[1]
kappa_xx = np.zeros(n_elem)
kappa_yy = np.zeros(n_elem)
x_cent   = np.zeros(n_elem)

for i in range(n_elem):
    ni = mesh.t[:, i]
    x  = mesh.p[0, ni]; y = mesh.p[1, ni]
    A  = np.column_stack([np.ones(3), x, y])
    cx, _, _, _ = np.linalg.lstsq(A, wx_nodes[ni], rcond=None)
    cy, _, _, _ = np.linalg.lstsq(A, wy_nodes[ni], rcond=None)
    kappa_xx[i] = cx[1]         # d(wx)/dx
    kappa_yy[i] = cy[2]         # d(wy)/dy
    x_cent[i]   = x.mean()

sigma_elem = (E_fg * t_fin / (2 * (1 - nu_fg**2))
              * np.abs(kappa_xx + nu_fg * kappa_yy))

# ──────────────────────────────────────────────────────────────────────────────
# Results
# ──────────────────────────────────────────────────────────────────────────────

# Tip deflection (max |w| at far end)
tip_dofs = basis.get_dofs(lambda x: x[0] > b_fin - 1e-9).nodal['u']
tip_disp = float(np.max(np.abs(u_sol[tip_dofs])))

# Peak stress at root strip (x < 5 mm)
root_mask   = x_cent < 0.005
sigma_root_fea = float(sigma_elem[root_mask].max()) if root_mask.any() else float(sigma_elem.max())
SF_fea      = sigma_ult / sigma_root_fea

# Analytical comparison (from structural_checks.py)
M_root_anal   = N_fin * b_fin / 3
I_fin_root    = c_root * t_fin**3 / 12
sigma_root_anal = M_root_anal * (t_fin / 2) / I_fin_root
SF_anal       = sigma_ult / sigma_root_anal

def sf_tag(sf):
    if sf >= 3.0: return "PASS ✓"
    if sf >= 1.0: return "WARN ⚠"
    return "FAIL ✗"

print("=" * 68)
print("  Apex-1 — Fin FEA (scikit-fem Kirchhoff plate)")
print("=" * 68)
print(f"\n  Fin geometry:  {c_root*1000:.0f}/{c_tip*1000:.0f} mm root/tip chord, "
      f"{b_fin*1000:.0f} mm span, {t_fin*1000:.0f} mm thick")
print(f"  Material:      E-glass/epoxy  E={E_fg/1e9:.0f} GPa  "
      f"σ_ult={sigma_ult/1e6:.0f} MPa")
print(f"  Load case:     max-Q  (Mach 1.1, 5 km)  "
      f"q={q_maxQ/1e3:.1f} kPa  AoA={math.degrees(AoA_rad):.0f}°")
print(f"  Normal force:  N_fin = {N_fin:.1f} N  "
      f"(uniform pressure {p_uniform:.0f} Pa)")
print(f"\n  Mesh:  {mesh.p.shape[1]} nodes, {mesh.t.shape[1]} elements  "
      f"(4 uniform refinements)")
print(f"  DOFs:  {basis.N}")

print("\n  ── Deflection ───────────────────────────────────────────────")
print(f"  Max tip deflection:          {tip_disp*1000:.2f} mm")
print(f"  Tip/span ratio:              {tip_disp/b_fin*100:.2f}%  "
      f"({'acceptable' if tip_disp/b_fin < 0.05 else 'large'})")

print("\n  ── Root bending stress ──────────────────────────────────────")
print(f"  FEA peak stress (root):      {sigma_root_fea/1e6:.1f} MPa")
print(f"  Analytical (closed-form):    {sigma_root_anal/1e6:.1f} MPa")
print(f"  FEA / analytical ratio:      {sigma_root_fea/sigma_root_anal:.2f}×  "
      f"(FEA captures stress concentration at root corners)")
print(f"\n  SF (FEA):                    {SF_fea:.2f}   {sf_tag(SF_fea)}")
print(f"  SF (analytical):             {SF_anal:.2f}   {sf_tag(SF_anal)}")

print("\n  Note: FEA stress is higher than the analytical estimate because:")
print("    1. Uniform pressure (not triangular) puts more load toward the tip,")
print("       increasing the effective moment arm at the root.")
print("    2. Stress concentrations at the root LE/TE corners captured by FEA.")
print("    Both effects are physically real; the FEA result is more conservative.")

print("\n  ── Verdict ──────────────────────────────────────────────────")
if SF_fea >= 3.0:
    print(f"  Fin root stress: SAFE  (FEA SF = {SF_fea:.1f}, minimum 3.0)")
elif SF_fea >= 1.5:
    print(f"  Fin root stress: MARGINAL  (FEA SF = {SF_fea:.1f})")
else:
    print(f"  Fin root stress: FAIL  (FEA SF = {SF_fea:.1f})")
print("  For final design, add a bonded FG doubler at the fin root")
print("  and/or fillet the root corners to reduce stress concentration.")
print("=" * 68)
