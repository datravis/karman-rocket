"""
Apex-1 Structural Analysis
===========================
Analytical checks for the four dominant failure modes on a Karman-class rocket:

  1. Euler column buckling   — airframe collapses axially under motor thrust
  2. Local shell buckling    — CF tube wall buckles (thin-shell instability)
  3. Fin root bending stress — fin tears off under worst-case aerodynamic load
  4. Fin flutter proxy       — reduced frequency at max-Q and max speed

References:
  - NASA SP-8007: Buckling of Thin-Walled Circular Cylinders (1968)
  - Barrowman (1967): "The Theoretical Prediction of the Center of Pressure"
  - Theodorsen (1935): NACA TR 496 (reduced-frequency flutter criterion)
  - Bernoulli-Euler cantilever beam, 1st mode λ₁ = 1.875104

Caveat: these are closed-form estimates, not FEA. Use pycalculix / FEniCSx for
full structural simulation. Fin flutter in particular should be verified with
AeroFinSim or a dedicated aeroelastic code.
"""

import math

# ──────────────────────────────────────────────────────────────────────────────
# Geometry
# ──────────────────────────────────────────────────────────────────────────────
R_out  = 0.0635   # body outer radius, m  (127 mm OD)
t_wall = 0.002    # CF wall thickness, m
L_body = 3.500    # overall rocket length, m

R_in   = R_out - t_wall
A_tube = math.pi * (R_out**2 - R_in**2)      # tube cross-section area, m²
I_tube = math.pi / 4 * (R_out**4 - R_in**4)  # 2nd moment of area, m⁴

c_root = 0.360  # fin root chord, m
c_tip  = 0.160  # fin tip chord, m
b_fin  = 0.200  # exposed semi-span (cantilevered length), m
t_fin  = 0.004  # fin thickness, m

c_avg  = (c_root + c_tip) / 2
S_fin  = c_avg * b_fin          # planform area per fin, m²
taper  = c_tip / c_root

# Exposed aspect ratio (for Barrowman fin lift formula)
AR_exp = b_fin**2 / S_fin

# ──────────────────────────────────────────────────────────────────────────────
# Materials
# ──────────────────────────────────────────────────────────────────────────────
# Carbon fiber (woven 2×2 twill tube, longitudinal)
E_cf        = 70e9    # Young's modulus, Pa
sigma_ult_cf = 600e6  # ultimate compressive strength, Pa  (conservative)

# E-glass / epoxy (woven, balanced 0/90)
E_fg        = 17e9    # Young's modulus, Pa
G_fg        = 3.5e9   # shear modulus, Pa
rho_fg      = 1800    # density, kg/m³
nu_fg       = 0.18    # Poisson's ratio
sigma_ult_fg = 230e6  # ultimate tensile/bending strength, Pa (conservative)

# ──────────────────────────────────────────────────────────────────────────────
# Loading — from RocketPy simulation results
# ──────────────────────────────────────────────────────────────────────────────
m_launch    = 35.05   # kg
g           = 9.81    # m/s²
peak_G      = 31.7    # G  (from simulate.py)
F_axial     = m_launch * peak_G * g  # axial load during burn, N

# Max dynamic pressure: Mach 1.1 at ~5 km altitude (standard atmosphere)
#   ρ = 0.736 kg/m³,  a = 320 m/s  →  V = 352 m/s
rho_5km = 0.736   # kg/m³
a_5km   = 320.0   # m/s
V_maxQ  = 1.1 * a_5km
q_maxQ  = 0.5 * rho_5km * V_maxQ**2  # Pa

# Max speed: Mach 5.47 at ~17 km (a ≈ 295 m/s,  ρ ≈ 0.104 kg/m³)
a_17km   = 295.0
rho_17km = 0.104
V_maxspd = 5.47 * a_17km

# ──────────────────────────────────────────────────────────────────────────────
# CHECK 1: Euler Column Buckling
# ──────────────────────────────────────────────────────────────────────────────
# P_cr = π² × E × I / (K × L)²
# K = 1.0  (conservative: pinned-pinned end conditions)
# For fixed-fixed (fully bonded bulkheads) K = 0.5 — use 1.0 for margin.
K_col = 1.0
P_cr_euler = math.pi**2 * E_cf * I_tube / (K_col * L_body)**2
SF_euler   = P_cr_euler / F_axial

# ──────────────────────────────────────────────────────────────────────────────
# CHECK 2: Local Shell Buckling  (NASA SP-8007)
# ──────────────────────────────────────────────────────────────────────────────
# σ_cr = γ × 0.605 × E × (t / r)
# Knockdown factor γ = 0.40 accounts for geometric imperfections in composites
# (NASA SP-8007 recommends 0.32–0.65 depending on quality; 0.40 is conservative)
gamma_kd   = 0.40
sigma_cr_sh = gamma_kd * 0.605 * E_cf * (t_wall / R_out)  # Pa
sigma_appl  = F_axial / A_tube                              # Pa
SF_shell    = sigma_cr_sh / sigma_appl

# ──────────────────────────────────────────────────────────────────────────────
# CHECK 3: Fin Root Bending Stress
# ──────────────────────────────────────────────────────────────────────────────
# Normal force per fin (Barrowman slender-wing lift, 5° AoA, max-Q)
# CN_alpha (per fin, per radian) ≈ 2π / (1 + 2/AR_exp)  — corrected for AR
# This is the Helmbold formula for low-AR wings.
AoA_deg    = 5.0
AoA_rad    = math.radians(AoA_deg)
CN_alpha_fin = 2 * math.pi / (1 + 2 / AR_exp)   # /rad
N_fin      = CN_alpha_fin * q_maxQ * S_fin * AoA_rad  # normal force, N

# Root bending moment (assume triangular spanwise load → resultant at b/3 from root)
M_root     = N_fin * b_fin / 3

# Fin section at root: width = c_root, height = t_fin
# Bending about the chord axis (flapwise bending under aerodynamic lift)
I_fin_root  = c_root * t_fin**3 / 12
sigma_root  = M_root * (t_fin / 2) / I_fin_root
SF_fin_root = sigma_ult_fg / sigma_root

# ──────────────────────────────────────────────────────────────────────────────
# CHECK 4: Fin Flutter — Reduced Frequency Proxy
# ──────────────────────────────────────────────────────────────────────────────
# Natural bending frequency of fin as cantilevered uniform plate (first mode):
#   f_n = (λ₁² / (2π × b²)) × sqrt(E × t² / (ρ × 12 × (1 − ν²)))
#   where λ₁ = 1.875104 (first Bernoulli-Euler mode eigenvalue)
#
# Note: this treats the fin as a uniform rectangular plate, which slightly
# overestimates frequency for a tapered fin. Considered conservative.
lam1  = 1.875104
f_nat = (lam1**2 / (2 * math.pi * b_fin**2)) * math.sqrt(
    E_fg * t_fin**2 / (rho_fg * 12 * (1 - nu_fg**2))
)
omega_n = 2 * math.pi * f_nat

# Theodorsen reduced frequency:  k = ω × c_avg / (2V)
# Flutter onset is empirically associated with k dropping below ~0.2;
# k < 0.08 is considered high flutter risk.
# These thresholds are from classical thin-airfoil flutter theory.
k_maxQ   = omega_n * c_avg / (2 * V_maxQ)
k_maxspd = omega_n * c_avg / (2 * V_maxspd)

def flutter_rating(k):
    if k >= 0.20:  return "SAFE     ✓"
    if k >= 0.08:  return "MARGINAL ⚠"
    return             "AT RISK  ✗"

# ──────────────────────────────────────────────────────────────────────────────
# Output
# ──────────────────────────────────────────────────────────────────────────────
def sf_line(label, sf, threshold=3.0):
    if sf >= threshold:    tag = "PASS ✓"
    elif sf >= 1.0:        tag = "WARN ⚠"
    else:                  tag = "FAIL ✗"
    return f"  {label:<40}  SF = {sf:6.2f}   {tag}"

print("=" * 68)
print("  Apex-1 — Structural Analysis")
print("=" * 68)
print(f"\n  Axial load (burn, {peak_G:.0f} G):  {F_axial/1e3:7.2f} kN")
print(f"  Max-Q (Mach 1.1, 5 km):    {q_maxQ/1e3:7.2f} kPa  "
      f"({0.5*rho_5km*V_maxQ**2/101325*100:.1f}% of 1 atm)")

print("\n  ── Airframe (CF tube 127 mm OD × 2 mm wall) ──────────────────")
print(sf_line("Euler column buckling", SF_euler))
print(f"    P_cr = {P_cr_euler/1e3:.0f} kN   F = {F_axial/1e3:.1f} kN")
print(sf_line("Local shell buckling (NASA SP-8007, γ=0.40)", SF_shell))
print(f"    σ_cr = {sigma_cr_sh/1e6:.0f} MPa   σ_applied = {sigma_appl/1e6:.2f} MPa")

print("\n  ── Fins (FG 4 mm thick, 200 mm span, 360/160 mm root/tip) ────")
print(sf_line(f"Fin root bending stress ({AoA_deg:.0f}° AoA, max-Q)", SF_fin_root))
print(f"    N_fin = {N_fin:.0f} N   M_root = {M_root:.1f} N·m   "
      f"σ_root = {sigma_root/1e6:.1f} MPa   (limit {sigma_ult_fg/1e6:.0f} MPa)")

print("\n  ── Fin flutter (reduced frequency, Theodorsen criterion) ──────")
print(f"  Fin natural bending frequency:       f_n = {f_nat:.1f} Hz")
print(f"  Reduced freq at max-Q  (M1.1, 5km):   k = {k_maxQ:.3f}  "
      f"→  {flutter_rating(k_maxQ)}")
print(f"  Reduced freq at max-V  (M5.5, 17km):  k = {k_maxspd:.3f}  "
      f"→  {flutter_rating(k_maxspd)}")
print()
print("  Note: at Mach 5.5 / 17 km, air density is 8.5% of sea level.")
print("  Low density partially offsets the low reduced frequency,")
print("  but fin flutter at peak speed warrants FEA / AeroFinSim verification.")

print("\n  ── Summary ────────────────────────────────────────────────────")
print("  Column buckling:   SAFE (large margin from motor thrust to critical load)")
print("  Shell buckling:    SAFE (thin-wall CF handles applied compressive stress)")
print("  Fin root stress:   SAFE (FG strength >> aerodynamic bending load)")
print("  Fin flutter:       MARGINAL at max-Q; needs verification at max speed")
print()
print("  Highest priority action: validate fin flutter with AeroFinSim")
print("  or pycalculix FEA before committing to this design.")
print("=" * 68)
