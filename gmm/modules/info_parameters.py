
# Dictionary containing parameters information
info_parameters = {
    # Centerlines parameters
    "node_diameters": {
        "units": "m",
        "units_conversion": 1,
        "3d_plot": {
            "title": "D nodes",
            "scale": 1e-6
        },
        "histogram": {
            "xlabel": "Capillary voxel diameter",
            "ylabel": "Number of capillary voxels"
        }
    },
    # Deposition parameters
    "deposition_number_clogged_inlet_links": {
        "units": "-",
        "units_conversion": 1,
        "2d_plot": {
            "title": "Clogged inlet links number",
            "scale": 1,
            "xlabel": "Time [s]",
            "ylabel": "Number of inlet links"
        }
    },
    "deposition_number_clogged_links": {
        "units": "-",
        "units_conversion": 1,
        "2d_plot": {
            "title": "Clogged links number",
            "scale": 1,
            "xlabel": "Time [s]",
            "ylabel": "Number of links"
        }
    },
    "deposition_onset": {
        "units": "-",
        "units_conversion": 1,
        "3d_plot": {
            "title": "Deposition",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Potential deposition",
            "ylabel": "Number of capillary voxels"
        }
    },
    "deposition_rate": {
        "units": "kg/m2",
        "units_conversion": 1,
        "3d_plot": {
            "title": "r_dep",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Deposition rate",
            "ylabel": "Number of capillary voxels"
        }
    },
    "deposition_shear_stress": {
        "units": "Pa",
        "units_conversion": 1,
        "3d_plot": {
            "title": "tau_dep",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Deposition shear stress",
            "ylabel": "Number of capillary voxels"
        }
    },
    # Erosion parameters
    "erosion_onset": {
        "units": "-",
        "units_conversion": 1,
        "3d_plot": {
            "title": "Erosion onset",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Potential erosion",
            "ylabel": "Number of capillary voxels"
        }
    },
    "erosion_rate": {
        "units": "kg/m2",
        "units_conversion": 1,
        "3d_plot": {
            "title": "r_er",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Erosion rate",
            "ylabel": "Number of capillary voxels"
        }
    },
    "erosion_shear_stress": {
        "units": "Pa",
        "units_conversion": 1,
        "3d_plot": {
            "title": "tau_er",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Erosion shear stress",
            "ylabel": "Number of capillary voxels"
        }
    },
    "erosion_time_scale": {
        "units": "s",
        "units_conversion": 1,
        "3d_plot": {
            "title": "t_er",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Erosion time scale",
            "ylabel": "Number of capillary voxels"
        }
    },

    # Flow parameters
    "maximum_flow_rate": {
        "units": "nL/s",
        "units_conversion": 1e12,
        "2d_plot": {
            "title": "Maximum flow rate evolution",
            "scale": 1,
            "xlabel": "Maximum flow rate",
            "ylabel": "Accumulated volume simulation [m^3]"
        }
    },
    "pressure_gradients": {
        "units": "Pa",
        "units_conversion": 1,
        "3d_plot": {
            "title": "Delta P",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Pressure gradients",
            "ylabel": "Number of capillary voxels"
        }
    },
    "reynolds_number": {
        "units": "-",
        "units_conversion": 1,
        "3d_plot": {
            "title": "Re",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Reynolds number",
            "ylabel": "Number of capillary voxels"
        }
    },
    "wall_shear_stress": {
        "units": "Pa",
        "units_conversion": 1,
        "3d_plot": {
            "title": "wall_shear_stress",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Wall shear stress",
            "ylabel": "Number of capillary voxels"
        }
    },

    # Geometry parameters
    "accumulated_volume_simulation": {
        "units": "nl",
        "units_conversion": 1e12,
        "2d_plot": {
            "title": "Accumulated (during the simulation) volume evolution",
            "scale": 1e12,
            "xlabel": "Time [s]",
            "ylabel": "Accumulated volume simulation [nl]"
        }
    },
    "accumulated_volume_time_step": {
        "units": "nl",
        "units_conversion": 1e12,
        "2d_plot": {
            "title": "Accumulated (during the time step) volume evolution",
            "scale": 1e12,
            "xlabel": "Time [s]",
            "ylabel": "Accumulated volume time step [nl]"
        }
    },
    "aspect_ratio": {
        "units": "-",
        "units_conversion": 1,
        "3d_plot": {
            "title": "Aspect ratio",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Aspect ratio",
            "ylabel": "Number of capillary voxels"
        }
    },
    "link_diameters": {
        "units": "m",
        "units_conversion": 1,
        "3d_plot": {
            "title": "D",
            "scale": 1e-6
        },
        "histogram": {
            "xlabel": "Capillary voxel diameter",
            "ylabel": "Number of capillary voxels"
        }
    },
    "link_length": {
        "units": "m",
        "units_conversion": 1,
        "3d_plot": {
            "title": "L",
            "scale": 1e-6
        },
        "histogram": {
            "xlabel": "Capillary voxel length",
            "ylabel": "Number of capillary voxels"
        }
    },
    "pore_volume": {
        "units": "[-]",
        "units_conversion": 1e12,
        "2d_plot": {
            "title": "Pore volume evolution",
            "scale": 1e12,
            "xlabel": "Time [s]",
            "ylabel": "Pore volume [nl]"
        }
    },
    "porosity": {
        "units": "[-]",
        "units_conversion": 1,
        "2d_plot": {
            "title": "Porosity evolution",
            "scale": 1,
            "xlabel": "Time [s]",
            "ylabel": "Porosity [-]"
        }
    },
    "reactive_area": {
        "units": "m^2",
        "units_conversion": 1,
        "3d_plot": {
            "title": "A_r",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Reactive area [m^2]",
            "ylabel": "Number of capillary voxels"
        }
    },
    "total_accumulated_volume_simulation": {
        "units": "nl",
        "units_conversion": 1e12,
        "2d_plot": {
            "title": "Total accumulated volume simulation evolution",
            "scale": 1e12,
            "xlabel": "Time [s]",
            "ylabel": "Total accumulated volume [nl]"
        }
    },
    "total_accumulated_volume_time_step": {
        "units": "nl",
        "units_conversion": 1e12,
        "2d_plot": {
            "title": "Total accumulated volume time step evolution",
            "scale": 1e12,
            "xlabel": "Time [s]",
            "ylabel": "Total accumulated volume [nl]"
        }
    },
    "void_space_volume": {
        "units": "nl",
        "units_conversion": 1e12,
        "2d_plot": {
            "title": "Void space volume evolution",
            "scale": 1e12,
            "xlabel": "Time [s]",
            "ylabel": "Void space volume [nl]"
        }
    },
    "void_space_volume_ratio": {
        "units": "-",
        "units_conversion": 1,
        "2d_plot": {
            "title": "phi(t)/phi(0) evolution",
            "scale": 1,
            "xlabel": "Time [s]",
            "ylabel": "Void space volume ratio [-]"
        }
    },

    # Precipitation parameters
    "maximum_precipitation_rate": {
        "units": "mol/(s m^2)",
        "units_conversion": 1,
        "2d_plot": {
            "title": "Max r_ppt evolution",
            "scale": 1,
            "xlabel": "Time [s]",
            "ylabel": "Maximum precipitation rate [mol/(s m^2)]"
        }
    },
    "precipitation_number_clogged_inlet_links": {
        "units": "-",
        "units_conversion": 1,
        "2d_plot": {
            "title": "Clogged inlet links number",
            "scale": 1,
            "xlabel": "Time [s]",
            "ylabel": "Number of inlet links"
        }
    },
    "precipitation_number_clogged_links": {
        "units": "-",
        "units_conversion": 1,
        "2d_plot": {
            "title": "Clogged links number",
            "scale": 1,
            "xlabel": "Time [s]",
            "ylabel": "Number of links"
        }
    },
    "precipitation_clogging_evolution": {
        "units": "-",
        "units_conversion": 1,
        "3d_plot": {
            "title": "Clogged links distribution",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Number of clogged links",
            "ylabel": "Number of capillary voxels"
        }
    },
    "precipitation_rate": {
        "units": "mol/(s m^2)",
        "units_conversion": 1,
        "3d_plot": {
            "title": "Precipitation rate",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Precipitation rate",
            "ylabel": "Number of capillary voxels"
        }
    },

    # Static parameters
    "flow_rate": {
        "units": "nL/s",
        "units_conversion": 1e12,
        "3d_plot": {
            "title": "Q",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Flow rate",
            "ylabel": "Number of capillary voxels"
        }
    },
    "flow_speed": {
        "units": "mm/s",
        "units_conversion": 1e3,
        "3d_plot": {
            "title": "V",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Flow speed",
            "ylabel": "Number of capillary voxels"
        }
    },
    "permeability": {
        "units": "-",
        "units_conversion": 1,
        "2d_plot": {
            "title": "Permeability",
            "scale": 1,
            "xlabel": "Time [s]",
            "ylabel": "Permeability [-]"
        }
    },
    "pressures": {
        "units": "Pa",
        "units_conversion": 1,
        "3d_plot": {
            "title": "P - min(P)",
            "scale": 1
        },
        "histogram": {
            "xlabel": "Pressures",
            "ylabel": "Number of capillary voxels"
        }
    }
}
