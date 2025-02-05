"""
Solar System Simulation with Tidal Graphs and Additional Information:
   - Displays the current Moon phase (New Moon, Waxing Crescent, First Quarter, etc.)
   - Indicates the current season (Spring, Summer, Autumn, Winter)
   - Gravitational force on Earth is now computed dynamically using the instantaneous
     distances to the Moon and the Sun.
   - (No world map pinpointing is performed.)
   
To run:
    pip install vpython
    python solar_sim_improved_gravity.py
"""

from vpython import *
import math

# ----------------------------------------------------
# Global Simulation Settings and Physical Constants
# ----------------------------------------------------
scale_mode = "useful"       # "useful" for visualization; "actual" for near–real scales
simulation_speed = 1.0      # simulation days advanced per loop iteration
paused = False
previous_speed = simulation_speed
dt = 0.1                    # base time step in simulation days per iteration

# Orbital periods and angular speeds (radians per day)
earth_orbit_omega    = 2 * math.pi / 365.0   # Earth's orbit around Sun (365 days)
earth_rotation_omega = 2 * math.pi / 1.0     # Earth's rotation (1 day)
moon_orbit_omega     = 2 * math.pi / 27.3    # Moon's orbit period (~27.3 days)

# Initial angles (in radians)
earth_angle = 0
moon_angle  = 0
# Track Earth's rotation (for geographic coordinate transformation)
earth_self_rotation = 0

# Physical constants
G = 6.67430e-11        # gravitational constant, m^3 kg^-1 s^-2
g = 9.81               # gravitational acceleration, m/s^2
mass_sun   = 1.989e30  # kg
mass_earth = 5.972e24  # kg
mass_moon  = 7.34767309e22  # kg

# Distances and radii (km and m)
R_earth_km = 149600000   # Earth-Sun average distance (km)
r_moon_km  = 384400      # Moon's average distance from Earth (km)
R_earth_m  = 6.371e6     # Earth's radius (m)

# For "actual" scale rendering, compress km → VPython units so that Earth's orbit is ~10 units.
actual_distance_scale = 10 / R_earth_km

# VPython scale parameters for two modes (sizes in VPython units)
scale_params = {
    "actual": {
        "earth_orbit": R_earth_km * actual_distance_scale,
        "sun_radius":    696340   * actual_distance_scale,
        "earth_radius":  6371     * actual_distance_scale,
        # (For the Moon we no longer use a constant orbit distance.)
        "moon_radius":   1737     * actual_distance_scale
    },
    "useful": {
        "earth_orbit": 10,
        "sun_radius":    2,
        "earth_radius":  0.5,
        # The "useful" average distance for the Moon (for scaling) is 2 units.
        "moon_orbit":    2,
        "moon_radius":   0.15
    }
}

# Parameters for the Moon’s elliptical orbit
a_moon = r_moon_km   # semi-major axis in km
e_moon = 0.0549      # eccentricity
i_moon = math.radians(5)  # inclination of 5°

# -------------------------------------
# Global Variables for Tidal Graphs and Data
# -------------------------------------
# Default coordinate (degrees) for tide plotting (user adjustable)
chosen_lon = 0.0  # longitude (degrees)
chosen_lat = 0.0  # latitude (degrees)

# Plot window for tide graph (in days)
plot_window = 7   # default: 7 days

# Flags to toggle tide curves
show_moon_tide     = True
show_sun_tide      = True
show_combined_tide = True

# Cumulative distances (km)
earth_distance_travelled = 0.0         # Earth relative to Sun
moon_distance_travelled_sun = 0.0      # Moon relative to Sun
moon_distance_travelled_earth = 0.0    # Moon relative to Earth

# For finite-difference estimation of Moon’s traveled distance and speed:
prev_moon_physical = None  # previous Moon position (in m, relative to Sun)
prev_moon_offset_m = None  # previous Moon offset (in m, relative to Earth)
prev_moon_offset_m_velocity = None  # for velocity estimate relative to Earth
prev_moon_physical_velocity = None    # for velocity estimate relative to Sun

# -------------------------------------
# Helper Functions for Tidal Calculations
# -------------------------------------
def tide_height_component(M, d, cos_theta):
    """
    Compute the equilibrium tide height (in m) from a tide-raising body.
      M: mass (kg)
      d: distance from Earth's center (m)
      cos_theta: cosine of the angle between local vertical and direction to the body.
    """
    return (G * M / (d**3)) * (R_earth_m**2) / g * (3 * cos_theta**2 - 1)/2

def compute_tide_at_point(point, earth_center, moon_pos, sun_pos):
    """
    Compute the equilibrium tide height (m) at a given point on Earth's surface.
      point: inertial position (m)
      earth_center: Earth's center (m)
      moon_pos: Moon's position (m)
      sun_pos: Sun's position (m)
    Returns a tuple: (h_moon, h_sun, h_total)
    """
    local_vertical = (point - earth_center).norm()
    # Moon tide
    cos_theta_moon = dot(local_vertical, (moon_pos - earth_center).norm())
    h_moon = tide_height_component(mass_moon, r_moon_km * 1e3, cos_theta_moon)
    # Sun tide
    cos_theta_sun = dot(local_vertical, (sun_pos - earth_center).norm())
    h_sun = tide_height_component(mass_sun, R_earth_km * 1e3, cos_theta_sun)
    return h_moon, h_sun, (h_moon + h_sun)

def get_earth_surface_point(lon_deg, lat_deg, earth_center, rotation_angle):
    """
    Compute the inertial position of a point on Earth's surface given its geographic
    longitude and latitude (in degrees) and the Earth’s rotation angle.
    """
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    r_local = R_earth_m * vector(math.cos(lat)*math.cos(lon),
                                 math.sin(lat),
                                 math.cos(lat)*math.sin(lon))
    r_inertial = rotate(r_local, angle=rotation_angle, axis=vector(0,1,0))
    return earth_center + r_inertial

# -------------------------------------
# Create the VPython Scene and Celestial Bodies
# -------------------------------------
scene = canvas(title="Solar System Simulation: Sun, Earth, Moon",
               width=800, height=600,
               center=vector(0,0,0),
               background=color.black)

sun = sphere(pos=vector(0,0,0),
             radius=scale_params[scale_mode]["sun_radius"],
             color=color.yellow,
             emissive=True)

earth = sphere(pos=vector(scale_params[scale_mode]["earth_orbit"], 0, 0),
               radius=scale_params[scale_mode]["earth_radius"],
               texture=textures.earth,
               make_trail=True,
               trail_type="curve",
               retain=100)

# (The Moon’s radius remains fixed; its position is updated below.)
moon = sphere(pos=earth.pos + vector(scale_params["useful"]["moon_orbit"], 0, 0),
              radius=scale_params[scale_mode]["moon_radius"],
              color=color.white,
              make_trail=True,
              trail_type="curve",
              retain=50)

# An arrow marker on Earth to visualize its rotation.
earth_axis = arrow(pos=earth.pos, 
                   axis=vector(earth.radius, 0, 0),
                   color=color.cyan)

# -------------------------------------
# GUI Callbacks for Solar System Controls
# -------------------------------------
def set_actual_scale(evt):
    global scale_mode
    scale_mode = "actual"
    sun.radius = scale_params["actual"]["sun_radius"]
    earth.pos = vector(scale_params["actual"]["earth_orbit"] * math.cos(earth_angle),
                       0,
                       scale_params["actual"]["earth_orbit"] * math.sin(earth_angle))
    earth.radius = scale_params["actual"]["earth_radius"]
    earth_axis.pos = earth.pos
    earth_axis.axis = vector(earth.radius, 0, 0)

def set_useful_scale(evt):
    global scale_mode
    scale_mode = "useful"
    sun.radius = scale_params["useful"]["sun_radius"]
    earth.pos = vector(scale_params["useful"]["earth_orbit"] * math.cos(earth_angle),
                       0,
                       scale_params["useful"]["earth_orbit"] * math.sin(earth_angle))
    earth.radius = scale_params["useful"]["earth_radius"]
    earth_axis.pos = earth.pos
    earth_axis.axis = vector(earth.radius, 0, 0)

def pause_resume(evt):
    global simulation_speed, paused, previous_speed
    if not paused:
        previous_speed = simulation_speed
        simulation_speed = 0
        paused = True
        pause_button.text = "Resume"
    else:
        simulation_speed = previous_speed
        paused = False
        pause_button.text = "Pause"

def speed_up(evt):
    global simulation_speed, paused, previous_speed
    if not paused:
        simulation_speed *= 2
        previous_speed = simulation_speed
    else:
        previous_speed *= 2
    scene.title = f"Solar System Simulation (Speed: {simulation_speed:.2f} days/loop)"

def slow_down(evt):
    global simulation_speed, paused, previous_speed
    if simulation_speed > 0.01:
        simulation_speed /= 2
        previous_speed = simulation_speed
    scene.title = f"Solar System Simulation (Speed: {simulation_speed:.2f} days/loop)"

def view_sun(evt):
    scene.center = sun.pos

def view_earth(evt):
    scene.center = earth.pos

# Create buttons for solar system controls
button_actual     = button(text="Actual Scale", bind=set_actual_scale)
button_useful     = button(text="Useful Scale", bind=set_useful_scale)
pause_button      = button(text="Pause", bind=pause_resume)
button_speed_up   = button(text="Speed Up", bind=speed_up)
button_slow_down  = button(text="Slow Down", bind=slow_down)
button_view_sun   = button(text="View Sun", bind=view_sun)
button_view_earth = button(text="View Earth", bind=view_earth)

# -------------------------------------
# Tidal Data and Graphing Interface
# -------------------------------------
info_text = wtext(text="")

wtext(text="\n\nEnter coordinate for tide graph (deg): ")
wtext(text="Longitude: ")
input_lon = winput(text="0", bind=lambda evt: None)
wtext(text=" Latitude: ")
input_lat = winput(text="0", bind=lambda evt: None)

def update_coordinates(evt):
    global chosen_lon, chosen_lat
    try:
        chosen_lon = float(input_lon.text)
        chosen_lat = float(input_lat.text)
    except:
        pass
button_update_coord = button(text="Update Coordinates", bind=update_coordinates)

# Buttons to change the moving window of the tide graph
def set_window_1_day(evt):
    global plot_window
    plot_window = 1

def set_window_1_week(evt):
    global plot_window
    plot_window = 7

def set_window_1_month(evt):
    global plot_window
    plot_window = 30

wtext(text="\nTime Window for Tide Graph: ")
button_window_day   = button(text="1 Day", bind=set_window_1_day)
button_window_week  = button(text="1 Week", bind=set_window_1_week)
button_window_month = button(text="1 Month", bind=set_window_1_month)

# Checkboxes to toggle tide curves
def toggle_moon_tide(c):
    global show_moon_tide
    show_moon_tide = c.checked

def toggle_sun_tide(c):
    global show_sun_tide
    show_sun_tide = c.checked

def toggle_combined_tide(c):
    global show_combined_tide
    show_combined_tide = c.checked

wtext(text="\nSelect Tide Curves to Display: ")
checkbox_moon = checkbox(text="Moon Tide", bind=toggle_moon_tide, checked=True)
checkbox_sun = checkbox(text="Sun Tide", bind=toggle_sun_tide, checked=True)
checkbox_combined = checkbox(text="Combined Tide", bind=toggle_combined_tide, checked=True)

# Create the graph for tide heights
tide_graph = graph(title="Tide Height vs Time", xtitle="Time (days)", ytitle="Tide Height (m)",
                   fast=False, width=600, height=300)
moon_tide_curve = gcurve(graph=tide_graph, color=color.blue, label="Moon Tide")
sun_tide_curve = gcurve(graph=tide_graph, color=color.red, label="Sun Tide")
combined_tide_curve = gcurve(graph=tide_graph, color=color.green, label="Combined Tide")

# -------------------------------------
# Main Simulation Loop
# -------------------------------------
t = 0  # simulation time in days

while True:
    rate(100)  # up to 100 iterations per second
    current_dt = simulation_speed * dt   # effective time step (in days)
    t += current_dt

    # Update orbital angles and Earth's rotation angle
    earth_angle += earth_orbit_omega * current_dt
    moon_angle  += moon_orbit_omega * current_dt
    earth_self_rotation += earth_rotation_omega * current_dt

    # Update Earth's position (for VPython rendering)
    if scale_mode == "actual":
        orbit_radius = scale_params["actual"]["earth_orbit"]
    else:
        orbit_radius = scale_params["useful"]["earth_orbit"]
    earth.pos = vector(orbit_radius * math.cos(earth_angle),
                       0,
                       orbit_radius * math.sin(earth_angle))
    # (Rotate Earth to simulate day-night cycle)
    earth.rotate(angle=earth_rotation_omega * current_dt, axis=vector(0,1,0))
    earth_axis.pos = earth.pos
    earth_axis.rotate(angle=earth_rotation_omega * current_dt, axis=vector(0,1,0))

    # -------------------------
    # Compute Physical Positions (in SI units)
    # -------------------------
    R_earth_m_si = R_earth_km * 1e3
    earth_physical = vector(R_earth_m_si * math.cos(earth_angle),
                            0,
                            R_earth_m_si * math.sin(earth_angle))
    sun_physical = vector(0, 0, 0)  # Sun is at the origin

    # ----------------------------------------
    # Compute Moon’s Position Using an Elliptical, Inclined Orbit
    # ----------------------------------------
    # Use the current moon_angle as an approximate true anomaly (ν)
    r_moon_current = a_moon * (1 - e_moon**2) / (1 + e_moon * math.cos(moon_angle))  # in km

    # In the Moon’s orbital plane (before inclination), the offset (in km) is:
    moon_offset_km = vector(r_moon_current * math.cos(moon_angle),
                            0,
                            r_moon_current * math.sin(moon_angle))
    # Rotate this offset about the x-axis by i_moon (5° inclination)
    moon_offset_km = rotate(moon_offset_km, angle=i_moon, axis=vector(1,0,0))

    # For VPython rendering:
    if scale_mode == "actual":
        moon_vpython_offset = moon_offset_km * actual_distance_scale
    else:
        # In "useful" mode, scale so that when r_moon_current equals the average, the offset is as before.
        scale_factor = scale_params["useful"]["moon_orbit"] / r_moon_km
        moon_vpython_offset = moon_offset_km * scale_factor
    moon.pos = earth.pos + moon_vpython_offset

    # For physical (SI) computations, convert the Moon offset from km to m:
    moon_offset_m = moon_offset_km * 1e3
    moon_physical = earth_physical + moon_offset_m

    # -------------------------
    # Finite-Difference Estimates for Cumulative Distance and Speeds
    # -------------------------
    # For Moon relative to Earth
    if prev_moon_offset_m is None:
        prev_moon_offset_m = moon_offset_m
    delta_moon_earth = mag(moon_offset_m - prev_moon_offset_m)
    moon_distance_travelled_earth += delta_moon_earth / 1000.0
    prev_moon_offset_m = moon_offset_m

    # Approximate instantaneous Moon speed (relative to Earth) using finite differences
    if prev_moon_offset_m_velocity is None:
        prev_moon_offset_m_velocity = moon_offset_m
    if current_dt * 86400 > 0:
        v_moon_rel_approx = (moon_offset_m - prev_moon_offset_m_velocity) / (current_dt * 86400)
    else:
        v_moon_rel_approx = vector(0,0,0)
    speed_moon_rel = mag(v_moon_rel_approx) / 1000.0  # km/s
    prev_moon_offset_m_velocity = moon_offset_m

    # For Moon speed relative to Sun
    if prev_moon_physical_velocity is None:
        prev_moon_physical_velocity = moon_physical
    if current_dt * 86400 > 0:
        v_moon_global_approx = (moon_physical - prev_moon_physical_velocity) / (current_dt * 86400)
    else:
        v_moon_global_approx = vector(0,0,0)
    speed_moon_global = mag(v_moon_global_approx) / 1000.0  # km/s
    prev_moon_physical_velocity = moon_physical

    # For Earth, we continue to use the analytical circular-orbit speed:
    v_earth = vector(-R_earth_m_si * math.sin(earth_angle),
                     0,
                     R_earth_m_si * math.cos(earth_angle)) * (earth_orbit_omega / 86400)
    speed_earth = mag(v_earth) / 1000.0  # km/s

    # -------------------------
    # Compute Global Tide Extremes via a Coarse Grid
    # -------------------------
    max_tide = -1e9
    min_tide = 1e9
    for lat_deg in range(-90, 91, 10):
        for lon_deg in range(0, 360, 10):
            pt = get_earth_surface_point(lon_deg, lat_deg, earth_physical, earth_self_rotation)
            h_m, h_s, h_tot = compute_tide_at_point(pt, earth_physical, moon_physical, sun_physical)
            max_tide = max(max_tide, h_tot)
            min_tide = min(min_tide, h_tot)

    # -------------------------
    # Compute Tide at the Chosen Geographic Coordinate
    # -------------------------
    chosen_point = get_earth_surface_point(chosen_lon, chosen_lat, earth_physical, earth_self_rotation)
    h_m_chosen, h_s_chosen, h_tot_chosen = compute_tide_at_point(chosen_point, earth_physical, moon_physical, sun_physical)

    # -------------------------
    # Compute Moon Phase (text label)
    # -------------------------
    # Compute the ecliptic longitudes (projected onto the x-z plane) for Moon and Sun as seen from Earth.
    lon_moon = math.degrees(math.atan2((moon_physical - earth_physical).z, (moon_physical - earth_physical).x))
    lon_sun = math.degrees(math.atan2((sun_physical - earth_physical).z, (sun_physical - earth_physical).x))
    raw_diff = lon_moon - lon_sun
    if raw_diff < 0:
        raw_diff += 360
    d_angle = raw_diff if raw_diff <= 180 else 360 - raw_diff

    if d_angle < 22.5:
        moon_phase = "New Moon"
    elif d_angle < 67.5:
        moon_phase = "Waxing Crescent" if raw_diff <= 180 else "Waning Crescent"
    elif d_angle < 112.5:
        moon_phase = "First Quarter" if raw_diff <= 180 else "Last Quarter"
    elif d_angle < 157.5:
        moon_phase = "Waxing Gibbous" if raw_diff <= 180 else "Waning Gibbous"
    else:
        moon_phase = "Full Moon"

    # -------------------------
    # Compute Season
    # -------------------------
    # Assume t=0 corresponds to the vernal equinox (approx. March 21).
    # An offset of 80 days is added so that (t + 80) mod 365 gives a day-of-year:
    #   0-93: Spring, 93-186: Summer, 186-279: Autumn, 279-365: Winter.
    offset = 80
    day_of_year = (t + offset) % 365
    if day_of_year < 93:
        season = "Spring"
    elif day_of_year < 186:
        season = "Summer"
    elif day_of_year < 279:
        season = "Autumn"
    else:
        season = "Winter"

    # -------------------------
    # Update the Tide Graph (with a moving x-axis window)
    # -------------------------
    tide_graph.xmin = t - plot_window
    tide_graph.xmax = t
    if show_moon_tide:
        moon_tide_curve.plot(t, h_m_chosen)
    if show_sun_tide:
        sun_tide_curve.plot(t, h_s_chosen)
    if show_combined_tide:
        combined_tide_curve.plot(t, h_tot_chosen)

    # -------------------------
    # Update the On-Screen Info Text
    # -------------------------
    # Compute dynamic gravitational forces using instantaneous distances.
    d_earth_sun = mag(earth_physical - sun_physical)
    F_sun = G * mass_sun * mass_earth / (d_earth_sun**2)
    d_earth_moon = mag(moon_physical - earth_physical)
    F_moon = G * mass_moon * mass_earth / (d_earth_moon**2)

    years_passed = t / 365.0
    info_text.text = f"""
Time: {t:.2f} days  ({years_passed:.2f} years)
Moon Phase: {moon_phase}
Season: {season}

Speeds (km/s):
    Earth relative to Sun: {speed_earth:.2f}
    Moon relative to Sun:  {speed_moon_global:.2f}
    Moon relative to Earth: {speed_moon_rel:.2f}

Gravitational Forces on Earth:
    From Sun: {F_sun:.2e} N
    From Moon: {F_moon:.2e} N

Coordinates (in km):
    Sun:   {vector(0,0,0)/1000}
    Earth: {earth_physical/1000}
    Moon:  {moon_physical/1000}

Distances (km):
    Sun-Earth: {mag(earth_physical)/1000:.2f}
    Earth-Moon: {mag(moon_physical - earth_physical)/1000:.2f}
    Moon-Sun:  {mag(moon_physical)/1000:.2f}

Tide Extremes (Equilibrium Tide Model):
    Maximum tide height over Earth: {max_tide:+.3f} m
    Minimum tide height over Earth: {min_tide:+.3f} m

At chosen coordinate (lon={chosen_lon:.1f}°, lat={chosen_lat:.1f}°):
    Moon tide: {h_m_chosen:+.3f} m
    Sun tide:  {h_s_chosen:+.3f} m
    Combined:  {h_tot_chosen:+.3f} m

Cumulative Distances (km):
    Earth traveled (relative to Sun): {earth_distance_travelled:,.0f}
    Moon traveled (relative to Earth): {moon_distance_travelled_earth:,.0f}

Assumptions:
  - Equilibrium tide on a spherical Earth.
  - Circular Earth orbit.
  - The Moon’s orbit is elliptical (e = 0.0549) and inclined by 5°.
  - Tidal contributions from Moon and Sun add linearly.
    (Dynamic effects, local geography, and resonances are neglected.)
    """

    # Update Earth's cumulative travel along its orbit (arc length)
    earth_distance_travelled += R_earth_km * (earth_orbit_omega * current_dt)
