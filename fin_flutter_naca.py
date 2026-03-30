"""
Apex-1 Fin Flutter Analysis — Piston-Theory Torsional Divergence
=================================================================
Implements the NACA TN 4197 piston-theory approach for fin flutter/divergence:

  • Integrates a 1D vertical ascent trajectory (same parameters as simulate.py)
    to obtain the altitude–speed–Mach–density profile through burnout.
  • At each point, computes the torsional divergence speed using Collar's
    cantilever divergence criterion with piston-theory lift slope:

      q_div = (π/2)² × GJ / (b² × CLα_strip × e × c̄)

    where:
      GJ        = G × c̄ × t³/3  (Saint-Venant torsional stiffness)
      CLα_strip = 4/M            (piston theory, supersonic M > 1.2)
                = 2π             (thin-airfoil, subsonic M < 0.8)
      e         = c̄/4           (shear-centre to aerodynamic-centre offset,
                                  flat-plate: SC at c̄/2, AC at c̄/4)

  • Reports results for two CLα models:
      (A) 2D strip (no AR correction) — conservative upper bound on loads
      (B) Helmbold AR-corrected      — better for low-AR fins (AR = 0.77)

  • Finds the altitude of minimum flutter safety factor.
  • Calculates the fin thickness required for SF_AR ≥ 1.5.

Note on 1D vs 6DOF:
  The 1D model gives peak speed ~1852 m/s; RocketPy 6DOF gives 1614 m/s.
  The 1D model is ~15% fast because it ignores trajectory curvature and
  non-vertical thrust components. Results for both speeds are reported.

References:
  - Collar (1946): aeroelastic triangle; torsional divergence criterion
  - NACA TN 4197 (1957): flutter of low-aspect-ratio lifting surfaces
  - Piston theory: CLα = 4/M for M > 1.2 (Lighthill, 1953)
  - Barrowman / Helmbold: finite-span CLα correction

Dependencies:
  pip install scipy
"""

import math
import numpy as np
from scipy.integrate import solve_ivp

# ──────────────────────────────────────────────────────────────────────────────
# Fin geometry and material  (matches structural_checks.py and simulate.py)
# ──────────────────────────────────────────────────────────────────────────────
b_fin  = 0.200    # exposed semi-span, m
c_root = 0.360    # root chord, m
c_tip  = 0.160    # tip chord, m
t_fin  = 0.004    # thickness, m

c_avg  = (c_root + c_tip) / 2   # mean chord, m
S_fin  = c_avg * b_fin           # planform area per fin, m²
AR_exp = b_fin**2 / S_fin        # exposed aspect ratio

G_fg   = 3.5e9    # shear modulus, E-glass/epoxy, Pa

# Saint-Venant torsional constant (thin flat plate, t << c):  J = c·t³/3
J_fin  = c_avg * t_fin**3 / 3
GJ     = G_fg * J_fin            # torsional stiffness, N·m²

# Aerodynamic centre offset from elastic axis (shear centre)
# Flat plate: AC at c/4, SC at c/2  →  e = c/4 (destabilising: AC forward of SC)
e_ac   = c_avg / 4

# AR correction denominator (Helmbold / low-AR subsonic analogy)
AR_factor = 1 + 2 / AR_exp

# ──────────────────────────────────────────────────────────────────────────────
# International Standard Atmosphere
# ──────────────────────────────────────────────────────────────────────────────
def isa(h_m):
    h = max(h_m, 0.0)
    if h < 11_000:
        T = 288.15 - 0.0065 * h
        p = 101_325 * (T / 288.15) ** 5.2561
    elif h < 20_000:
        T = 216.65
        p = 22_632 * math.exp(-9.81 * (h - 11_000) / (287 * T))
    else:
        T = 216.65 + 0.001 * (h - 20_000)
        p = 5_474.9 * math.exp(-9.81 * (h - 20_000) / (287 * T))
    rho = p / (287 * T)
    a_s = math.sqrt(1.4 * 287 * T)
    return rho, p, a_s

# ──────────────────────────────────────────────────────────────────────────────
# 1D ascent trajectory  (same parameters as simulate.py)
# ──────────────────────────────────────────────────────────────────────────────
ELEV    = 1191.0   # m ASL  (Black Rock Desert launch site)
m_launch = 35.05   # kg
m_prop   = 23.75   # kg propellant
t_burn   = 10.0    # s
thrust   = 5400.0  # N (constant — neutral BATES burn)
A_body   = math.pi * 0.0635**2   # reference area, m²

# Mach-dependent drag coefficient (from simulate.py drag_profile)
_DRAG_PTS = [(0,.38),(0.7,.38),(0.85,.42),(0.95,.55),(1.05,.65),
             (1.15,.62),(1.30,.55),(1.60,.48),(2.00,.40),(2.50,.34),
             (3.00,.30),(4.00,.25),(5.00,.21),(7.00,.17),(10.0,.14)]

def cd_mach(M):
    if M <= _DRAG_PTS[0][0]:  return _DRAG_PTS[0][1]
    if M >= _DRAG_PTS[-1][0]: return _DRAG_PTS[-1][1]
    for i in range(len(_DRAG_PTS) - 1):
        m0, m1 = _DRAG_PTS[i][0], _DRAG_PTS[i+1][0]
        if m0 <= M <= m1:
            f = (M - m0) / (m1 - m0)
            return _DRAG_PTS[i][1] + f * (_DRAG_PTS[i+1][1] - _DRAG_PTS[i][1])

def odes(t, y):
    h, v = y
    rho, _, a_s = isa(h)
    mass  = m_launch - min(t, t_burn) * (m_prop / t_burn)
    thr   = thrust if t <= t_burn else 0.0
    q_dyn = 0.5 * rho * v * abs(v)
    drag  = cd_mach(abs(v) / a_s) * q_dyn * A_body
    return [v, (thr - drag) / mass - 9.81]

print("Integrating 1D ascent trajectory …", flush=True)
sol = solve_ivp(odes, [0, 300], [ELEV, 0.1], max_step=0.05, dense_output=True)

# Sample ascent up to apogee
t_apo_idx = int(np.argmax(sol.y[0]))
t_end     = sol.t[t_apo_idx]
t_pts     = np.linspace(0.5, t_end, 5000)
h_pts     = np.array([sol.sol(t)[0] for t in t_pts])
v_pts     = np.array([sol.sol(t)[1] for t in t_pts])

atm_pts   = [isa(h) for h in h_pts]
rho_pts   = np.array([x[0] for x in atm_pts])
a_s_pts   = np.array([x[2] for x in atm_pts])
M_pts     = v_pts / a_s_pts
q_dyn_pts = 0.5 * rho_pts * v_pts**2

# Key trajectory statistics
idx_maxv  = int(np.argmax(v_pts))
idx_maxq  = int(np.argmax(q_dyn_pts))

# ──────────────────────────────────────────────────────────────────────────────
# Section lift-curve slope
# ──────────────────────────────────────────────────────────────────────────────
def cla_section(M, ar_correct=False):
    """2D section CLα.  Blended subsonic→supersonic across transonic."""
    if M > 1.10:
        cla = 4.0 / M                        # piston theory
    elif M > 0.85:
        cla_sub = 2 * math.pi
        cla_sup = 4.0 / 1.10
        cla = cla_sub + (cla_sup - cla_sub) * (M - 0.85) / 0.25
    else:
        cla = 2 * math.pi                    # thin-airfoil, subsonic
    return cla / AR_factor if ar_correct else cla

# ──────────────────────────────────────────────────────────────────────────────
# Torsional divergence speed
# ──────────────────────────────────────────────────────────────────────────────
def v_diverge(M, rho, ar_correct=False):
    """Critical divergence speed (m/s) via Collar cantilever criterion."""
    cla  = cla_section(M, ar_correct)
    if cla <= 0:
        return float('inf')
    # q_div = (π/2)² · GJ / (b² · CLα · e · c̄)
    q_div = (math.pi / 2)**2 * GJ / (b_fin**2 * cla * e_ac * c_avg)
    return math.sqrt(2 * q_div / rho)

# ──────────────────────────────────────────────────────────────────────────────
# Sweep: flutter margin along ascent
# ──────────────────────────────────────────────────────────────────────────────
vd_2d = np.array([v_diverge(M, rho) for M, rho in zip(M_pts, rho_pts)])
vd_ar = np.array([v_diverge(M, rho, ar_correct=True) for M, rho in zip(M_pts, rho_pts)])

# Only consider supersonic flight (subsonic loads are lower; flutter more relevant ≥ M1.2)
super_mask = M_pts > 1.2

SF_2d = np.where(super_mask, vd_2d / v_pts, np.inf)
SF_ar = np.where(super_mask, vd_ar / v_pts, np.inf)

idx_min_ar = int(np.argmin(SF_ar))
idx_min_2d = int(np.argmin(SF_2d))

# ──────────────────────────────────────────────────────────────────────────────
# Required thickness for target safety factor
# ──────────────────────────────────────────────────────────────────────────────
def required_thickness(SF_target, v_rocket, M_rocket, rho_rocket, ar_correct=False):
    """Solve V_div = SF_target × V_rocket for t_fin.
       V_div ∝ t^(3/2)  →  t_req = t_fin × (SF_target × V_rocket / V_div_current)^(2/3)
    """
    cla = cla_section(M_rocket, ar_correct)
    # q_div_req = 0.5 × rho × (SF_target × V_rocket)²
    q_div_req = 0.5 * rho_rocket * (SF_target * v_rocket)**2
    # q_div = (π/2)² × G × c̄ × t³ / (3 × b² × cla × e × c̄)
    #       = (π/2)² × G × t³ / (3 × b² × cla × e)
    # Solve for t:
    t3 = q_div_req * 3 * b_fin**2 * cla * e_ac / ((math.pi/2)**2 * G_fg)
    return t3 ** (1/3)

# Critical condition: minimum AR-corrected SF
M_crit  = M_pts[idx_min_ar]
v_crit  = v_pts[idx_min_ar]
rho_crit = rho_pts[idx_min_ar]
h_crit  = h_pts[idx_min_ar]

t_req_15_ar = required_thickness(1.5, v_crit, M_crit, rho_crit, ar_correct=True)
t_req_30_ar = required_thickness(3.0, v_crit, M_crit, rho_crit, ar_correct=True)
t_req_15_2d = required_thickness(1.5, v_crit, M_crit, rho_crit, ar_correct=False)

# Also check against the known 6DOF peak speed (1614 m/s) at the same altitude
V_6DOF_peak = 1614.0   # m/s, from simulate.py
M_6DOF_peak = V_6DOF_peak / a_s_pts[idx_maxv]   # approx same altitude
SF_ar_6dof  = v_diverge(M_pts[idx_min_ar], rho_pts[idx_min_ar], ar_correct=True) / V_6DOF_peak

# ──────────────────────────────────────────────────────────────────────────────
# Output
# ──────────────────────────────────────────────────────────────────────────────
def sf_tag(sf):
    if sf >= 3.0: return "SAFE ✓"
    if sf >= 1.5: return "MARGINAL ⚠"
    if sf >= 1.0: return "AT RISK ✗"
    return             "DIVERGENCE ✗✗"

print("=" * 70)
print("  Apex-1 — Fin Flutter / Torsional Divergence  (NACA TN 4197)")
print("=" * 70)
print(f"\n  Fin:      {c_root*1000:.0f}/{c_tip*1000:.0f} mm root/tip chord  ·  "
      f"{b_fin*1000:.0f} mm span  ·  {t_fin*1000:.0f} mm thick")
print(f"  Material: E-glass/epoxy  G = {G_fg/1e9:.1f} GPa")
print(f"  GJ:       {GJ:.4f} N·m²  (Saint-Venant, thin plate)")
print(f"  AR_exp:   {AR_exp:.3f}   AR-correction factor: {AR_factor:.2f}×")

print("\n  ── Trajectory  (1D vertical, same parameters as simulate.py) ────")
print(f"  Apogee (1D):     {h_pts[-1]/1000:.1f} km ASL  "
      f"(RocketPy 6DOF: 109.5 km)")
print(f"  Peak speed (1D): {v_pts[idx_maxv]:.0f} m/s  Mach {M_pts[idx_maxv]:.2f}  "
      f"at {h_pts[idx_maxv]/1000:.1f} km  ← burnout")
print(f"  Peak speed (6DOF, from simulate.py): {V_6DOF_peak:.0f} m/s  Mach 5.47")
print(f"  Max-Q (1D):      {q_dyn_pts[idx_maxq]/1e3:.0f} kPa  "
      f"at {h_pts[idx_maxq]/1000:.1f} km  Mach {M_pts[idx_maxq]:.2f}")
print(f"  Note: 1D peak speed ~15% higher than 6DOF (ignores trajectory curvature)")

print("\n  ── Divergence speed vs rocket speed  (supersonic phase) ──────────")
print(f"  {'Alt km':>7}  {'M':>5}  {'V rkt':>7}  {'Vd 2D':>8}  {'SF 2D':>7}  "
      f"{'Vd AR':>8}  {'SF AR':>7}")
print("  " + "-"*65)

# Print ~15 representative points in the supersonic phase
sup_idxs = np.where(super_mask)[0]
step = max(1, len(sup_idxs) // 12)
shown = list(sup_idxs[::step]) + [idx_min_ar]
shown = sorted(set(shown))
for i in shown:
    tag = " ← min SF" if i == idx_min_ar else ""
    print(f"  {h_pts[i]/1000:>7.1f}  {M_pts[i]:>5.2f}  "
          f"{v_pts[i]:>7.0f}  {vd_2d[i]:>8.0f}  {SF_2d[i]:>6.2f}   "
          f"{vd_ar[i]:>8.0f}  {SF_ar[i]:>6.2f}{tag}")

print("\n  ── Critical condition (minimum AR-corrected SF) ───────────────────")
print(f"  Altitude:     {h_crit/1000:.1f} km  (t ≈ {t_pts[idx_min_ar]:.1f} s)")
print(f"  Mach:         {M_crit:.2f}")
print(f"  V_rocket (1D): {v_crit:.0f} m/s     V_div_AR: {vd_ar[idx_min_ar]:.0f} m/s")
print(f"  SF (2D strip, conservative):  {SF_2d[idx_min_ar]:.2f}   {sf_tag(SF_2d[idx_min_ar])}")
print(f"  SF (AR-corrected, realistic): {SF_ar[idx_min_ar]:.2f}   {sf_tag(SF_ar[idx_min_ar])}")
print(f"\n  Adjusted for 6DOF peak speed ({V_6DOF_peak:.0f} m/s at same altitude):")
print(f"  SF (AR-corrected, 6DOF speed): {SF_ar_6dof:.2f}   {sf_tag(SF_ar_6dof)}")

print("\n  ── Required fin thickness ─────────────────────────────────────────")
print(f"  Current:           {t_fin*1000:.1f} mm  →  SF_AR = {SF_ar[idx_min_ar]:.2f}")
print(f"  For SF_AR ≥ 1.5:   {t_req_15_ar*1000:.1f} mm")
print(f"  For SF_AR ≥ 3.0:   {t_req_30_ar*1000:.1f} mm")
print(f"  (Conservative 2D, SF ≥ 1.5):  {t_req_15_2d*1000:.1f} mm")

print("\n  ── Physical interpretation ────────────────────────────────────────")
print("  Torsional divergence: aerodynamic pitching moment (about the fin's")
print("  shear centre at mid-chord) overcomes the fin's torsional stiffness.")
print("  For a symmetric flat-plate fin, this is the primary aeroelastic")
print("  instability — classical bending-torsion flutter requires CG/SC offset.")
print()
print("  The AR correction (Helmbold: CLα_eff = CLα_2D / (1 + 2/AR)) accounts")
print("  for the strong 3D tip relief on this low-AR fin (AR = 0.77).")
print("  At AR = 0.77, tip vortices reduce effective loads by ÷3.6×,")
print("  which is physically significant but should be verified by 3D CFD.")
print()
if SF_ar[idx_min_ar] < 1.0:
    print("  ⚠ Torsional divergence is predicted below peak rocket speed.")
    print("  Increasing fin thickness to ≥ {:.0f} mm is strongly recommended".format(t_req_15_ar*1000))
    print("  before committing to this design.")
elif SF_ar[idx_min_ar] < 1.5:
    print("  ⚠ Safety factor is MARGINAL. Verify with 3D CFD / FEA.")
    print(f"  Increasing fin thickness to {t_req_15_ar*1000:.0f} mm would give SF ≥ 1.5.")
else:
    print("  Fin is flutter-safe per this analysis.")
    print("  Recommend 3D CFD or AeroFinSim verification before final design lock.")

print("=" * 70)
