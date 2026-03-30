"""
Export Apex-1 rocket design to OpenRocket (.ork) and RASP (.eng) formats.

Outputs:
  Q54000.eng   — RASP thrust curve for the custom Q-class motor
  Apex1.ork    — OpenRocket file (import via File > Open)

To use:
  1. In OpenRocket: Edit > Preferences > User-defined thrust curves → Add Q54000.eng
  2. File > Open → Apex1.ork
  3. The motor "Q54000" will be selectable from the motor database
  4. Run the pre-configured "Black Rock Desert" simulation

Mass note: OpenRocket calculates component masses from geometry + density.
  Override them to match the RocketPy model if needed:
    Rocket dry (no motor):  5.80 kg,  CG 1.20 m from nose
    Motor hardware:         5.50 kg
    APCP propellant:       23.75 kg
    Launch mass:           35.05 kg
"""

import os
import uuid
import zipfile

# ── RASP Motor File ────────────────────────────────────────────────────────────
# Columns: name  diam_mm  len_mm  delays  prop_mass_kg  total_mass_kg  mfr
# Thrust pairs: time_s  thrust_N  (terminated by ';')
ENG_CONTENT = """\
; Apex-1 Q54000 — custom single-use APCP motor
; 6-grain BATES design, 118 mm OD casing, 1940 mm length
; Total impulse: 54,000 N·s  (Q class: 40,961–81,920 N·s)
; Average thrust: 5,400 N    Burn time: 10 s    Isp ≈ 232 s
; Propellant: 23.75 kg       Motor dry: 5.50 kg
;
Q54000 118 1940 0 23.75 29.25 Custom
   0.000     0.0
   0.100  5400.0
  10.000  5400.0
  10.100     0.0
;
"""

# ── OpenRocket XML ─────────────────────────────────────────────────────────────
# All linear dimensions in metres.
#
# Rocket layout (nose tip = 0):
#   0.000 – 0.700  Tangent ogive nose cone (CF)
#   0.700 – 3.500  Body tube 127 mm OD, 2 mm CF wall
#     1.560 – 3.500  Motor mount inner tube 118 mm OD (flush aft)
#     2.950          Fin root LE (4 × trapezoidal FG)
#                    root 360 mm | tip 160 mm | span 200 mm | sweep 130 mm
#   3.500            Nozzle exit

def _id() -> str:
    return str(uuid.uuid4())


def build_ork_xml() -> str:
    rid  = _id()   # rocket
    sid  = _id()   # stage
    nid  = _id()   # nose
    bid  = _id()   # body tube
    mid  = _id()   # motor mount inner tube
    fid  = _id()   # fin set
    did  = _id()   # drogue
    p1id = _id()   # main chute 1
    p2id = _id()   # main chute 2
    p3id = _id()   # main chute 3
    cfg  = _id()   # motor configuration

    # Body tube length = total length − nose cone = 3.500 − 0.700 = 2.800 m
    # Fin root LE from body tube front = 2.950 − 0.700 = 2.250 m
    # Motor mount (1.940 m) placed flush with aft end of body tube

    return f"""\
<?xml version='1.0' encoding='utf-8'?>
<openrocket version="1.7" creator="Apex-1 export script">

  <rocket>
    <name>Apex-1</name>
    <id>{rid}</id>
    <designer>Dan Travis</designer>
    <comment>Karman line attempt — 54,000 N·s Q-class APCP, 127 mm CF airframe.
IMPORTANT: Load Q54000.eng as a user-defined motor before simulating
(Edit › Preferences › User-defined thrust curves).</comment>

    <motorconfiguration configid="{cfg}" default="true">
      <name>Q54000 sustainer</name>
      <motor component="{mid}" configid="{cfg}">
        <designation>Q54000</designation>
        <ignitionevent>launch</ignitionevent>
        <ignitiondelay>0.0</ignitiondelay>
        <overhang>0.000</overhang>
      </motor>
    </motorconfiguration>

    <subcomponents>
      <stage>
        <name>Sustainer</name>
        <id>{sid}</id>
        <subcomponents>

          <!-- Nose cone: 700 mm tangent ogive, carbon fiber ─────── -->
          <nosecone>
            <name>Tangent Ogive Nose</name>
            <id>{nid}</id>
            <finish>normal</finish>
            <material type="bulk" density="1600.0">Carbon fiber</material>
            <length>0.700</length>
            <shape>ogive</shape>
            <shapeclipped>false</shapeclipped>
            <shapeparameter>1.0</shapeparameter>
            <aftradius>0.0635</aftradius>
            <aftshoulderradius>0.0615</aftshoulderradius>
            <aftshoulderlength>0.060</aftshoulderlength>
            <aftshoulderthickness>0.002</aftshoulderthickness>
            <aftshouldercapped>false</aftshouldercapped>
            <isflipped>false</isflipped>
          </nosecone>

          <!-- Body tube: 2800 mm, 127 mm OD, 2 mm CF wall ────────── -->
          <bodytube>
            <name>Airframe</name>
            <id>{bid}</id>
            <finish>normal</finish>
            <material type="bulk" density="1600.0">Carbon fiber</material>
            <length>2.800</length>
            <thickness>0.002</thickness>
            <radius>0.0635</radius>
            <subcomponents>

              <!-- Motor mount inner tube: 118 mm OD, flush aft ───── -->
              <innertube>
                <name>Motor Mount</name>
                <id>{mid}</id>
                <material type="bulk" density="1600.0">Carbon fiber</material>
                <length>1.940</length>
                <radialposition>0.0</radialposition>
                <radialdirection>0.0</radialdirection>
                <outerradius>0.059</outerradius>
                <thickness>0.003</thickness>
                <axialoffset method="bottom">0.000</axialoffset>
                <motormount>true</motormount>
              </innertube>

              <!-- 4 × trapezoidal FG fins ──────────────────────────
                   Root LE: 2.950 m from nose = 2.250 m from body tube front
                   Root chord: 360 mm | Tip chord: 160 mm
                   Span: 200 mm      | Sweep (LE): 130 mm          -->
              <trapezoidfinset>
                <name>Fins</name>
                <id>{fid}</id>
                <finish>normal</finish>
                <material type="bulk" density="1800.0">Fiberglass</material>
                <fincount>4</fincount>
                <rotation>0.0</rotation>
                <thickness>0.004</thickness>
                <crosssection>square</crosssection>
                <cant>0.0</cant>
                <rootchord>0.360</rootchord>
                <tipchord>0.160</tipchord>
                <height>0.200</height>
                <sweeplength>0.130</sweeplength>
                <axialoffset method="top">2.250</axialoffset>
              </trapezoidfinset>

              <!-- Drogue: 24" round, deploys at apogee ────────────
                   Cd·S = 0.75 × π/4 × 0.6096² = 0.275 m² ≈ 0.28 m²  -->
              <parachute>
                <name>Drogue 24in</name>
                <id>{did}</id>
                <finish>normal</finish>
                <material type="surface" density="0.060">Nylon</material>
                <cd>0.750</cd>
                <diameter>0.6096</diameter>
                <linecount>8</linecount>
                <linelength>0.900</linelength>
                <linematerial type="line" density="0.0018">Nylon</linematerial>
                <axialoffset method="top">0.100</axialoffset>
                <deployevent>apogee</deployevent>
                <deployaltitude>0.0</deployaltitude>
                <deploydelay>0.2</deploydelay>
              </parachute>

              <!-- Main cluster: 3 × 60" rounds, deploy at 300 m AGL
                   Each: Cd·S = 0.75 × π/4 × 1.524² = 1.368 m²
                   Total: 3 × 1.368 = 4.10 m²                       -->
              <parachute>
                <name>Main 1 of 3 (60in)</name>
                <id>{p1id}</id>
                <finish>normal</finish>
                <material type="surface" density="0.060">Nylon</material>
                <cd>0.750</cd>
                <diameter>1.524</diameter>
                <linecount>12</linecount>
                <linelength>3.000</linelength>
                <linematerial type="line" density="0.0018">Nylon</linematerial>
                <axialoffset method="top">0.200</axialoffset>
                <deployevent>altitude</deployevent>
                <deployaltitude>300.0</deployaltitude>
                <deploydelay>0.2</deploydelay>
              </parachute>

              <parachute>
                <name>Main 2 of 3 (60in)</name>
                <id>{p2id}</id>
                <finish>normal</finish>
                <material type="surface" density="0.060">Nylon</material>
                <cd>0.750</cd>
                <diameter>1.524</diameter>
                <linecount>12</linecount>
                <linelength>3.000</linelength>
                <linematerial type="line" density="0.0018">Nylon</linematerial>
                <axialoffset method="top">0.200</axialoffset>
                <deployevent>altitude</deployevent>
                <deployaltitude>300.0</deployaltitude>
                <deploydelay>0.2</deploydelay>
              </parachute>

              <parachute>
                <name>Main 3 of 3 (60in)</name>
                <id>{p3id}</id>
                <finish>normal</finish>
                <material type="surface" density="0.060">Nylon</material>
                <cd>0.750</cd>
                <diameter>1.524</diameter>
                <linecount>12</linecount>
                <linelength>3.000</linelength>
                <linematerial type="line" density="0.0018">Nylon</linematerial>
                <axialoffset method="top">0.200</axialoffset>
                <deployevent>altitude</deployevent>
                <deployaltitude>300.0</deployaltitude>
                <deploydelay>0.2</deploydelay>
              </parachute>

            </subcomponents>
          </bodytube>

        </subcomponents>
      </stage>
    </subcomponents>
  </rocket>

  <simulations>
    <simulation status="notrun">
      <name>Black Rock Desert — BALLS site</name>
      <simulator>RK4Simulator</simulator>
      <calculator>BarrowmanCalculator</calculator>
      <conditions>
        <configid>{cfg}</configid>
        <launchrodlength>6.0</launchrodlength>
        <launchrodangle>2.0</launchrodangle>
        <launchroddirection>0.0</launchroddirection>
        <windaverage>0.0</windaverage>
        <windturbulence>0.1</windturbulence>
        <launchaltitude>1191.0</launchaltitude>
        <launchlatitude>40.859</launchlatitude>
        <launchlongitude>-119.065</launchlongitude>
        <atmosphere model="isa"/>
        <timestep>0.05</timestep>
      </conditions>
    </simulation>
  </simulations>

</openrocket>
"""


def main() -> None:
    out_dir = os.path.dirname(os.path.abspath(__file__))

    # RASP .eng motor file
    eng_path = os.path.join(out_dir, "Q54000.eng")
    with open(eng_path, "w") as f:
        f.write(ENG_CONTENT)
    print(f"Motor file:      {eng_path}")

    # OpenRocket .ork file  (ZIP archive containing the XML)
    ork_path = os.path.join(out_dir, "Apex1.ork")
    with zipfile.ZipFile(ork_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("rocket.ork", build_ork_xml())
    print(f"OpenRocket file: {ork_path}")

    print()
    print("Usage:")
    print("  1. OpenRocket → Edit › Preferences › User-defined thrust curves")
    print("       Add: Q54000.eng")
    print("  2. File › Open → Apex1.ork")
    print("  3. Motor 'Q54000' will be available in the motor chooser")
    print("  4. Optionally override component masses to match RocketPy model:")
    print("       Rocket dry (no motor):  5.80 kg  (CG 1.20 m from nose)")
    print("       Motor hardware:         5.50 kg")
    print("       APCP propellant:       23.75 kg")
    print("       Launch mass:           35.05 kg")


if __name__ == "__main__":
    main()
