Solar System Simulation with Tidal Graphs
==========================================

Description:
------------
This project is an interactive simulation of a simplified solar system implemented in Python using VPython.
It simulates the following:

  - The Sun (displayed as a stationary, glowing yellow sphere)
  - The Earth (orbiting the Sun on a circular path while rotating on its axis)
  - The Moon (orbiting the Earth in an elliptical, 5°‑inclined orbit)

In addition, the simulation computes tidal heights based on a simplified equilibrium tide model.
Real‑time tide graphs are generated, displaying:
  • The tidal effect due to the Moon,
  • The tidal effect due to the Sun, and
  • Their combined effect,
all calculated at a user‑specified geographic coordinate.

The simulation also provides real‑time dynamic information including:
  • Simulation time (in days and years)
  • Speeds of the Earth and Moon
  • Instantaneous gravitational forces acting on the Earth
  • Cumulative distances traveled by the Earth and Moon
  • Current moon phase (e.g., New Moon, Full Moon, etc.)
  • Current season (Spring, Summer, Autumn, or Winter)

Features:
---------
• Celestial Bodies:
    – The Sun, Earth, and Moon are rendered in 3D.
    – The Earth’s rotation simulates the day–night cycle.
    – The Moon follows an elliptical orbit with eccentricity e = 0.0549 and is inclined by 5°.

• Tidal Calculations:
    – Tidal heights are computed using the formula:
      
         h = (G * M / d³) * (Rₑ² / g) * ((3 * cos²θ - 1) / 2)
      
      where:
        • G = gravitational constant,
        • M = mass of the tide‑raising body (Moon or Sun),
        • d = instantaneous distance from the Earth’s center,
        • Rₑ = Earth’s radius,
        • g = gravitational acceleration,
        • θ = angle between the local vertical and the direction to the tide‑raising body.

• Graphical Output:
    – A real‑time graph plots tide heights over a selectable time window.
    – On‑screen information updates dynamically (speeds, gravitational forces, distances, etc.).

• Interactive Controls:
    – Buttons allow the user to change simulation speed, pause/resume the simulation,
      and switch between different scaling modes.
    – Input fields let the user set the geographic coordinate for which tidal data is computed.

Assumptions:
------------
1. The Earth’s orbit around the Sun is circular.
2. The Earth is modeled as a perfect sphere.
3. The Moon’s orbit is elliptical (with e = 0.0549) and inclined by 5° relative to Earth’s orbital plane.
4. The tidal model assumes an instantaneous (equilibrium) response.
5. Tidal contributions from the Moon and Sun are added linearly, neglecting dynamic and local effects.

Installation and Usage:
-----------------------
1. Requirements:
   - Python 3.x
   - VPython (Install via pip: `pip install vpython`)

2. To Run the Simulation:
   - Clone this repository.
   - Open a terminal/command prompt in the project directory.
   - Execute the simulation script with:
     
         python solar_sim.py
     
   - A VPython window will open displaying the interactive simulation along with tidal graphs and controls.

Repository Structure:
---------------------
  - solar_sim.py      : Main Python simulation script.
  - README.txt        : This documentation file.
  - (Additional directories may include resources such as images, CSS, or JS files if integrated with a website.)

Notes:
------
- This simulation is designed for educational and demonstrative purposes.
- Some simplifications (e.g., equilibrium tide model, approximated moon phase calculations) have been made.
- Future improvements may include more accurate orbital dynamics, dynamic ocean responses, and further interactive features.

Author:
-------
Koen Ceton
