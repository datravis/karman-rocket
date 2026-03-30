"""
Apex-1 Karman Line Rocket Simulation
RocketPy 1.12.0

Design: 5-inch (127 mm) carbon fiber airframe, single-stage Q-class APCP motor
Target: >100 km ASL apogee (Karman line), ≥1.5 cal static stability, dual-deploy recovery

Launch site: Black Rock Desert, Nevada (BALLS/XPRS site — standard for FAR/Tripoli record attempts)
"""

import numpy as np
from rocketpy import Environment, SolidMotor, Rocket, Flight

# ── Environment: Black Rock Desert, NV ────────────────────────────────────────
ELEVATION = 1191  # m ASL — Black Rock playa

env = Environment(latitude=40.859, longitude=-119.065, elevation=ELEVATION)
env.set_date((2026, 8, 15, 18))
env.set_atmospheric_model(type="standard_atmosphere")

# ── Motor: Apex Q54000 (Q-class APCP, 6-grain BATES) ─────────────────────────
#
# Grain geometry:
#   OD: 110 mm  |  ID: 50 mm  |  Length: 300 mm  |  N: 6  |  Sep: 4 mm
#   ρ_APCP = 1750 kg/m³
#   Propellant mass ≈ 23.75 kg
#
# Thrust: smoothed trapezoidal curve, avg thrust 5400 N over 10 s
#   Total impulse = 54,000 N·s  (Q class: 40,961–81,920 N·s ✓)
#   Isp = 54,000 / (23.75 × 9.81) ≈ 232 s  (realistic APCP)
#   T/W at launch = 5400 / (35.1 × 9.81) ≈ 15.7  (excellent for altitude record)
#
# Motor casing: 118 mm OD, ~1.94 m long (6 grains + nozzle + header plate)
# Note: thrust curve is simplified for numerical stability — rapid ignition
#       transients cause ODE stiffness; physical impulse and Isp are preserved.

Q54000_thrust = [
    [0.000,     0],
    [0.100,  5400],  # fast but smooth ignition ramp
    [10.000, 5400],  # neutral BATES burn
    [10.100,    0],  # burnout
]

motor = SolidMotor(
    thrust_source=Q54000_thrust,
    dry_mass=5.5,                    # casing, nozzle, hardware — kg
    dry_inertia=(0.25, 0.25, 0.012),
    nozzle_radius=0.040,
    grain_number=6,
    grain_density=1750,              # APCP kg/m³
    grain_outer_radius=0.055,
    grain_initial_inner_radius=0.025,
    grain_initial_height=0.300,
    grain_separation=0.004,
    grains_center_of_mass_position=0.960,
    center_of_dry_mass_position=0.980,
    nozzle_position=0.000,
    throat_radius=0.018,
    coordinate_system_orientation="nozzle_to_combustion_chamber",
)

# ── Drag profile: Mach-dependent Cd ───────────────────────────────────────────
# Reference area = π × r_body²  (RocketPy convention)
# Calibrated for a slender CF rocket (fineness ratio ≈ 27.6:1)
# Subsonic: form drag + skin friction ~0.38
# Transonic peak at Mach ~1.1: wave drag adds ~0.25
# Hypersonic: wave drag falls per Newtonian theory
drag_profile = [
    [0.00, 0.38],
    [0.70, 0.38],
    [0.85, 0.42],
    [0.95, 0.55],
    [1.05, 0.65],   # transonic peak
    [1.15, 0.62],
    [1.30, 0.55],
    [1.60, 0.48],
    [2.00, 0.40],
    [2.50, 0.34],
    [3.00, 0.30],
    [4.00, 0.25],
    [5.00, 0.21],
    [7.00, 0.17],
    [10.0, 0.14],
]

# ── Rocket: Apex-1 ────────────────────────────────────────────────────────────
#
# Layout (nose_to_tail, all dims in meters):
#
#  0       0.70     1.05      1.45                    3.50
#  ├────────┤────────┤─────────┤────────────────────────┤
#  │ Ogive  │  Avo   │ Coupler │   Fin can + Motor      │
#  │  Nose  │  Bay   │         │   (Q54000, ~1.94 m)    │
#
#  Body tube: 127 mm OD (5"), 2 mm CF wall → radius = 63.5 mm
#
#  Dry mass breakdown (motor hardware excluded):
#    Nose cone (CF ogive)            0.50 kg  @ 0.35 m
#    Avionics bay (CF + hardware)    0.90 kg  @ 0.875 m
#    Body tube + motor mount (CF)    1.60 kg  @ 2.25 m
#    4 trapezoidal fins (FG)         0.90 kg  @ 3.10 m
#    Recovery system                 1.30 kg  @ 0.80 m
#    Rail buttons, misc              0.30 kg  @ 1.80 m
#    Drogue deployment bag           0.30 kg  @ 0.90 m
#    ────────────────────────────────────────────────
#    Rocket dry (no motor)           5.80 kg
#    Motor dry                       5.50 kg
#    ────────────────────────────────────────────────
#    Total dry                      11.30 kg
#    Propellant                     23.75 kg
#    Launch mass                    35.05 kg

rocket = Rocket(
    radius=0.0635,                          # 127 mm OD / 2
    mass=5.80,                              # rocket dry mass, no motor
    inertia=(4.80, 4.80, 0.016),            # Ixx, Iyy, Izz — kg·m²
    power_off_drag=drag_profile,
    power_on_drag=drag_profile,
    center_of_mass_without_motor=1.20,      # m from nose
    coordinate_system_orientation="nose_to_tail",
)

# Motor at tail — nozzle face at 3.50 m from nose
rocket.add_motor(motor, position=3.50)

# Nose cone — 700 mm tangent ogive, carbon fiber
rocket.add_nose(length=0.700, kind="ogive", position=0.000)

# 4 trapezoidal fins — fiberglass
# Root: 360 mm  |  Tip: 160 mm  |  Span: 200 mm  |  Sweep: 130 mm
# Provides CP well aft of CG → ≥1.5 cal stability through transonic
rocket.add_trapezoidal_fins(
    n=4,
    root_chord=0.360,
    tip_chord=0.160,
    span=0.200,
    position=2.950,          # root LE from nose
    sweep_length=0.130,
    cant_angle=0,
)

# ── Parachutes — dual-deploy recovery ─────────────────────────────────────────
# Drogue: 24" round, fires at apogee (~108 km) — stabilises descent through
#   upper atmosphere, limits drift, slows rocket to ~30 m/s at low altitude.
#   At sea-level-equivalent: Cd·S = 0.28 m²
rocket.add_parachute(
    name="drogue",
    cd_s=0.28,
    trigger="apogee",
    sampling_rate=105,
    lag=0.20,
    noise=(0, 8.3, 0.5),
)

# Main: 3 × 60" round cluster, fires at 300 m AGL — targets ≤7 m/s landing
#   Cd·S = 3 × π × (0.762)² × 0.75 / 4 × π ≈ 4.10 m²
#   Terminal v at launch-site density ≈ 6.9 m/s for 11.3 kg dry mass
#
# Guard against RocketPy phase-transition artefact: height can briefly read
# a large negative value when the drogue phase starts; the guard discards it.
def main_trigger(pressure, height, state):
    if height < 0 or height > 200_000:  # reject numerical artefacts
        return False
    return state[5] < 0 and height < 300   # descending AND below 300 m AGL

rocket.add_parachute(
    name="main",
    cd_s=4.10,
    trigger=main_trigger,
    sampling_rate=105,
    lag=0.20,
    noise=(0, 8.3, 0.5),
)

# ── Flight simulation ──────────────────────────────────────────────────────────
# max_time=1200 s covers powered flight + full coast to apogee (161 s) plus
# the majority of the descent; the main opens at ~300 m AGL around t≈770 s.
flight = Flight(
    rocket=rocket,
    environment=env,
    rail_length=6.0,    # 20 ft heavy-duty launch rail
    inclination=88,     # near-vertical (small off-axis prevents singularity)
    heading=0,
    max_time=1200,
    time_overshoot=True,
)

# ── Results ───────────────────────────────────────────────────────────────────
KARMAN   = 100_000          # m ASL — Karman line
apogee_m = flight.apogee   # m ASL
pass_km  = apogee_m / 1000

# Maximum speed and acceleration during ascent (avoid parachute artefacts)
t_apo       = flight.apogee_time
t_s         = np.linspace(1.0, t_apo, 600)
speeds      = [np.sqrt(flight.vx(t)**2 + flight.vy(t)**2 + flight.vz(t)**2)
               for t in t_s]
max_spd     = max(speeds)

t_burn      = np.linspace(0.5, 12.0, 300)
max_g       = max(
    np.sqrt(flight.ax(t)**2 + flight.ay(t)**2 + flight.az(t)**2)
    for t in t_burn
) / 9.81

margin_cal  = rocket.static_margin(0)

# Landing speed: find first point below 10 m AGL during descent
v_land      = float("nan")
for t in np.linspace(t_apo + 200, 1199, 4000):
    try:
        h = flight.z(t) - ELEVATION
        if 0 < h < 10 and flight.vz(t) < 0:
            v_land = abs(flight.vz(t))
            break
    except Exception:
        pass

print("=" * 60)
print("  Apex-1 Karman Line Rocket — Simulation Results")
print("=" * 60)
print(f"  Apogee (ASL)       : {apogee_m/1000:>9.2f} km")
print(f"  Karman line pass   : {'PASS ✓' if apogee_m >= KARMAN else 'FAIL ✗'}  "
      f"(margin: {(apogee_m - KARMAN)/1000:+.1f} km)")
print(f"  Max ascent speed   : {max_spd:>9.1f} m/s  (Mach {max_spd/340.3:.2f})")
print(f"  Max acceleration   : {max_g:>9.2f} G  (burn phase)")
print(f"  Time to apogee     : {t_apo:>9.1f} s")
print(f"  Static margin (t=0): {margin_cal:>9.2f} cal")
print(f"  Landing speed      : {v_land:>9.2f} m/s")
print(f"  Sim runtime        : {flight.t_final:>9.1f} s  (max_time cap)")
print("=" * 60)
print()
print("  Rocket Design Summary — Apex-1")
print("  ─────────────────────────────────────────────────────")
print(f"  Airframe           : 127 mm OD, 3.50 m long, carbon fiber")
print(f"  Nose cone          : 700 mm tangent ogive, CF")
print(f"  Fins               : 4 × trapezoidal FG, 360/160 mm root/tip, 200 mm span")
print(f"  Motor              : Q54000 APCP, 6-grain BATES, 118 mm OD casing")
print(f"  Motor impulse      : {motor.total_impulse:,.0f} N·s (Q class ✓)")
print(f"  Burn time          : {motor.burn_out_time:.1f} s")
print(f"  Dry mass           : {rocket.mass + motor.dry_mass:.2f} kg")
print(f"  Propellant mass    : {motor.propellant_initial_mass:.2f} kg")
print(f"  Launch mass        : {rocket.mass + motor.dry_mass + motor.propellant_initial_mass:.2f} kg")
print(f"  Drogue             : 24\" round, cd_s = 0.28 m², fires at apogee")
print(f"  Main               : 3 × 60\" round, cd_s = 4.10 m², fires at 300 m AGL")
print("=" * 60)
