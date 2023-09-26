Geometry Modification Module
----------------------------
The Geometry Modification Module (GMM) is a simulation toolkit for analyzing the geometry evolution of porous media due to pore-scale processes, considering the effects of flow conditions, geometry, and mineral reactions within a capillary network (See the [Box folder containing the prior art references, tutorials, and PowerPoint presentations](https://ibm.box.com/s/w2t6lbqjak3506zc67ud06jc9mtxmutm). The geometry evolution models and test case parameters are based on phenomenological correlations from the literature (See the [GMM Box](https://ibm.box.com/s/w2t6lbqjak3506zc67ud06jc9mtxmutm)). It is possible to add more sets of test cases or other correlations.

The main features of the GMM are:

1) A set of Python programs for simulating:
    - Spatiotemporal evolution of capillary networks due to pore-scale processes,
    - Mineral dissolution of calcite,
    - Mineral precipitation of calcite.
    - Clogging due to the accumulation of precipitated material within the capillary network.

2) A set of JSON files with test cases for simulating the geometry evolution of calcite rock samples (i.e., Berea_100x100x100, Berea_500x500x500, and MRS_100x100x400) and the effect of the precipitation rate on capillary clogging, including: 
    - Erosion according to the Erosion and Deposition Law in [Jager, R., et al. (2017)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.013110) and [Matias, A., et al. (2021)](https://doi.org/10.1016/j.jocs.2021.101360);
  
    - Mineral dissolution correlations in [Molins, S., et al. (2014)](https://pubs.acs.org/doi/10.1021/es5013438), [(2017)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR020323), [(2021)](https://link.springer.com/article/10.1007/s10596-019-09903-x), and
  
    - Mineral precipitation correlations of [Noiriel, C., et al. (2012)](https://www.sciencedirect.com/science/article/pii/S0009254112002331) and [Lasaga (1981)](https://www.degruyter.com/document/doi/10.1515/9781501508233-008/html). The following methods for simulating mineral precipitation within capillary networks:
      - Uniform mineral precipition: the precipitation rate remains constant during the simulations. The input parameters include: 
        - `Constant_precipitation_rate_input`: precipitation rate
        - `Constant_precipitation_rate_parameters`: a set of paratmters in [Noiriel, C., et al. (2012)](https://www.sciencedirect.com/science/article/pii/S0009254112002331) and [Lasaga (1981)](https://www.degruyter.com/document/doi/10.1515/9781501508233-008/html), including:
          - calcium_concentration_variation;
          - experiment_factor;
          - m_coefficient;
          - n_coefficient;
          - precipitation_rate_constant,
          - precipitation_reaction_time;
          - saturation_index

      - Flow dependent precipitation: the precipitation rate is a function of the flow-rate, calcium concentration variation between the inlet and the outlet of the rock sample, and the reactive area (see Eq. 3 in [Noiriel, C., et al. (2012)](https://www.sciencedirect.com/science/article/pii/S0009254112002331)).
        - `Constant_flow_rate`: a flow-dependent precipitation rate is calculated at the first tieration and kept constant;
        - `Flow_dependent_precipitation_rate`: a flow-dependent precipitation rate is calculated each time step.

3) A set of Python programs for creating:
    - A set of `centerlines.json` files containing:
       - Modified capillary network geometry: link diameter and link length distribution;
    - A set of `.h5` files containing: 
       - Erosion-related parameters:
         - condition for the existence of erosion,
         - erosion rate,
         - erosion shear stress, and
         - erosion time-scale;
  
      - Precipitation-related parameters:
         - precipitation number first clogged,
         - precipitation clogging evolution,
         - precipitation clogging onset, and
         - precipitation clogging time step;
  
      - Flow parameters:
        - pressure gradients,
        - reynolds number, and
        - wall shear stress;
  
      - Geometry parameters:
        - accumulated volume,
        - aspect ratio,
        - link length,
        - link radius,
        - porosity,
        - reactive area,
        - total accumulated volume,
        - void space volume,
        - void space volume ratio, and
        - void space volume variation;
  
       - Static results obtained from the Flow Discovery Simulator:
         - flow speed,
         - flow rate,
         - pressures, and 
         - permeability.

4) A Python program for creating 2d plots, 3d plots, and histograms of the parameters in (1):
    - For specific simulation time steps and
    - Ratios between two chosen time steps.
  
5) A Python program for creating `.mp4` videos and snapshots of the parameters in (1):
    - For specific simulation time steps and
    - Ratios between two chosen time steps.

6) An additional Python program in the `util/` folder for creating 3d plots and histograms of parameters related to the Flow Discovery Simulator (i.e., flow speed, flow rate, node diameters, link diameters, pressures, and pressure gradients).

Input parameters
----------------

The input parameters of the Model consist of a set of phasic (i.e., liquid and solid) properties, flow initial and boundary conditions, spatial domain characteristics, and pore-scale process parameters defined in a `.json` file. Also, the digital rock sample `centerlines.json` file or the `centerlines.json` file must be provided. The following input parameters may vary depending on the phenomenological correlation or geometry evolution model.

**Phasic, flow, and geometry**
| Parameter                | Units      | Type            |
|--------------------------|------------|-----------------|
| Liquid density           | kg m^-3    | Phasic property |
| Liquid dynamic viscosity | Pa s       | Phasic property |
| Solid density            | kg m^-3    | Phasic property |
| Solid molar mass         | kg mol^-1  | Phasic property |
| Porosity                 | -          | Geometry        |
| Capillary length         | voxel      | Geometry        |
| Capillary radius         | voxel      | Geometry        |
| Voxel size               | m          | Geometry        |
| Pressure                 | Pa         | Flow            |
| Flow rate                | m^3 s^-1   | Flow            |
| Flow speed               | m s^-1     | Flow            |
| Minimal diameter size    | m          | Simulation      |
| Reaction time            | s          | Simulation      |
| Simulation time          | s          | Simulation      |

**Pore-scale processes**
| Parameter                             | SI units          | Pore-scale process    |
|-------------------------------------  |-----------------  |--------------------   |
| Erodibility coefficient               | s m^-1            | Erosion               |
| Deposition coefficient                | s m^-1            | Deposition            |
| Relative concentration                | -                 | Deposition            |
| Inlet concentration                   | mol m^-3          | Dissolution           |
| Reaction constant                     | mol m^-1 s^-1     | Dissolution           |
| Thermodynamic activity of H^+         | m^-3 mol^-1       | Dissolution           |
| Calcite solubility                    |                   | Precipitation         |
| Thermodynamic activity of Ca^2+       | m^-3 mol^-1       | Precipitation         |
| Thermodynamic activity of CO~3^2+    | m^-3 mol^-1       | Precipitation         |

Output dataset parameters
-----------------

The outputs of the Model consist of a set of `.h5` files containing the following centerlines, flow, geometry, static, and pore-scale process parameters:
**Centerlines**
| Parameter                     | Units         |
|-----------------------------  |------------   |
| Updated link squared radius   | voxel^2       |

**Flow**
| Parameter                     | Units         |
|-----------------------------  |------------   |
| Pressure gradients            | Pa            |
| Reynolds number               | -             |
| Wall shear stress             | Pa            |

**Geometry**
| Parameter                     | Units         |
|-----------------------------  |------------   |
| Aspect ratio                  | -             |
| Porosity                      | -             |
| Reactive area                 | m^2           |
| Volume variation              | m^3           |
| Void space volume             | m^3           |
| Void space ratio              | -             |

**Static**
| Parameter                     | Units         |
|-----------------------------  |------------   |
| Flow rate                     | m^3 s^-1      |
| Flow speed                    | m s^-1        |
| Pressure                      | Pa            |

**Pore-scale processes**

| Parameter                         | Units             |
|-------------------------------    |-----------------  |
| Dissolution rate                  | mol m^-2 s^-1     |
| Erosion onset                     | -                 |
| Erosion rate                      | kg m^-2           |
| Erosion shear stress              | Pa                |
| Erosion time scale                | s                 |
| Clogged capillaries               | -                 |
| Clogging evolution                | -                 |
| Number first clogged capillaries  | -                 |
| Precipitation rate                | mol m^-2 s^-1     |
| Time for first clogging           | s                 |
