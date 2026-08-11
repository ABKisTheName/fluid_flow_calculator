"""
===============================================================================
ENGINEERING COMPUTATION MODULE
===============================================================================
Contains hydraulic calculations for incompressible pipe flow using:
  - Darcy-Weisbach Equation for pressure loss.
  - Swamee-Jain Equation for explicit friction factor estimation.
===============================================================================
"""

import math
import numpy as np


def convert_to_si(flow_rate_m3h, pipe_diameter_mm, roughness_mm, fluid_viscosity_cp):
    """Converts standard engineering input units to strict SI base units."""
    Q = flow_rate_m3h / 3600.0        # m³/s
    D = pipe_diameter_mm / 1000.0      # m
    e = roughness_mm / 1000.0          # m
    mu = fluid_viscosity_cp / 1000.0   # Pa·s
    return Q, D, e, mu


def calculate_flow_hydraulics(flow_rate_m3h, pipe_diameter_mm, roughness_mm, 
                              fluid_density, fluid_viscosity_cp, pipe_length_m):
    """
    Computes hydraulic velocity, Reynolds number, friction factor, 
    and pressure drop across a given pipe length.
    """
    # Unit Conversions
    Q, D, e, mu = convert_to_si(flow_rate_m3h, pipe_diameter_mm, roughness_mm, fluid_viscosity_cp)
    rho = fluid_density
    L = pipe_length_m

    # Area and Velocity
    area = (math.pi / 4.0) * (D ** 2)
    velocity = Q / area

    # Reynolds Number
    Re = (rho * velocity * D) / mu

    # Friction Factor & Regime Classification
    if Re < 2300:
        flow_regime = "Laminar Flow"
        f_d = 64.0 / Re
    elif 2300 <= Re <= 4000:
        flow_regime = "Transitional Flow"
        f_laminar = 64.0 / 2300
        f_turb = 0.25 / (math.log10((e / (3.7 * D)) + (5.74 / (4000 ** 0.9)))**2)  
        f_d = f_laminar + (Re - 2300) * (f_turb - f_laminar) / (4000 - 2300)
    else:
        flow_regime = "Turbulent Flow"
        f_d = 0.25 / (math.log10((e / (3.7 * D)) + (5.74 / (Re ** 0.9)))**2 )  

    # Darcy-Weisbach Pressure Drop: ΔP = f_d * (L/D) * (ρ * v² / 2)
    delta_p_pa = f_d * (L / D) * (rho * (velocity ** 2) / 2.0)
    delta_p_bar = delta_p_pa / 100000.0
    delta_p_psi = delta_p_pa / 6894.76

    return {
        "velocity": velocity,
        "reynolds_number": Re,
        "friction_factor": f_d,
        "flow_regime": flow_regime,
        "pressure_drop_pa": delta_p_pa,
        "pressure_drop_bar": delta_p_bar,
        "pressure_drop_psi": delta_p_psi
    }


def generate_diameter_sensitivity(flow_rate_m3h, current_diameter_mm, roughness_mm, 
                                   fluid_density, fluid_viscosity_cp, pipe_length_m, num_points=50):
    """Generates sensitivity data mapping pipe diameter variations to pressure drop."""
    min_d = max(10.0, current_diameter_mm * 0.25)
    max_d = current_diameter_mm * 2.0
    diameters_mm_range = np.linspace(min_d, max_d, num_points)
    dp_bar_list = []

    for d_mm in diameters_mm_range:
        res = calculate_flow_hydraulics(
            flow_rate_m3h, d_mm, roughness_mm, fluid_density, fluid_viscosity_cp, pipe_length_m
        )
        dp_bar_list.append(res["pressure_drop_bar"])

    return diameters_mm_range, dp_bar_list