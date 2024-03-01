FlowDiscovery FlowSimulator
=======================

[![Build Status](https://app.travis-ci.com/IBM/flowdiscovery-simulator.svg?branch=master)](https://app.travis-ci.com/IBM/flowdiscovery-simulator)

For developers
--------------

Read [CONTRIBUTING.md](CONTRIBUTING.md).

Command-line tools
------------------

| Name                                                    | Version      | License      |
|:-------------------------------------------------------:|:------------:|:------------:|
| [CMake](https://cmake.org)                              | `3.20.5`     | BSD-3-Clause |
| [Coverage.py](https://coverage.readthedocs.io)          | `5.5.0`      | Apache-2.0   |
| [Cpplint](https://github.com/cpplint/cpplint)           | `1.5.5`      | BSD-3-Clause |
| [Doxygen](http://www.doxygen.org)                       | `1.9.1`      | GPL-2.0      |
| [GCC](https://gcc.gnu.org)                              | `11.1.1`     | GPL-3.0      |
| [gperftools](https://github.com/gperftools/gperftools)  | `2.9.1`      | BSD-3-Clause |
| [LCOV](http://ltp.sourceforge.net/coverage)             | `1.14`       | GPL-2.0      |
| [Ninja](https://ninja-build.org)                        | `1.10.2`     | Apache-2.0   |
| [Pipenv](https://pipenv.pypa.io)                        | `2020.11.15` | MIT          |

C++ libraries
-------------

| Name                                                | Version   | License      |
|:---------------------------------------------------:|:---------:|:------------:|
| [Armadillo](http://arma.sourceforge.net)            | `10.6.0`  | Apache-2.0   |
| [GoogleLog](https://github.com/google/glog)         | `0.3.5`   | BSD-3-Clause |
| [GoogleTest](https://github.com/google/googletest)  | `1.10.0`  | BSD-3-Clause |
| [RapidJSON](http://rapidjson.org)                   | `1.1.0`   | MIT          |
| [SUNDIALS](https://github.com/LLNL/sundials)        | `5.6.1`   | BSD-3-Clause |
| [TCLAP](http://tclap.sourceforge.net)               | `1.2.3`   | MIT          |

Python libraries
----------------

| Name                                                                      | Version  | License      |
|:-------------------------------------------------------------------------:|:--------:|:------------:|
| [Celery](http://www.celeryproject.org)                                    | `5.1.2`  | BSD-3-Clause |
| [ibm-cos-sdk](https://github.com/ibm/ibm-cos-sdk-python/)                 | `2.10.0` | Apache-2.0   |
| [jsonschema](https://github.com/Julian/jsonschema)                        | `3.2.0`  | MIT          |
| [matplotlib](http://matplotlib.org)                                       | `3.4.2`  | Python-2.0   |
| [numpy-stl](https://github.com/WoLpH/numpy-stl)                           | `2.16.0` | BSD-3-Clause |
| [numpy](http://www.numpy.org)                                             | `1.21.1` | BSD-3-Clause |
| [OpenPNM](http://openpnm.org)                                             | `2.7.0`  | MIT          |
| [plumbum](https://plumbum.readthedocs.io)                                 | `1.7.0`  | MIT          |
| [scikit-image](http://scikit-image.org)                                   | `0.18.2` | BSD-3-Clause |
| [scipy](https://www.scipy.org)                                            | `1.7.0`  | BSD-3-Clause |

JSON Schemas
------------

| Name                                                | Version   | License      |
|:---------------------------------------------------:|:---------:|:------------:|
| [JSON Graph Format](http://jsongraphformat.info/)   | `3`       | MIT          |

License
------------

<!-- All source files must include a Copyright and License header. The SPDX license header is
preferred because it can be easily scanned. -->

If you would like to see the detailed LICENSE click [here](LICENSE).

```text
#
# Copyright 2020- IBM Inc. All rights reserved
# SPDX-License-Identifier: Apache-2.0
#
```

Related datasets
------------

**Digital Rock datasets**: Microtomography datasets available in the [Digital Rocks Portal](https://dx.doi.org/10.17612/f4h1-w124) and [Figshare](https://doi.org/10.25452/figshare.plus.21375565).

Related repositories
------------

**Digital Rock**: The code used to obtain the CNM representations of the rock samples is available at <https://github.com/IBM/flowdiscovery-digital-rock>. Additional algorithms used for processing and segmenting the raw grayscale images, are available as Python code at <https://github.com/IBM/microCT-Dataset>.
 
**Flow Simulator**: Code for running multi-phase flow simulations and the geometry evolution due to pore-scale processes in porous media modeled as network of capillaries is available at <https://github.com/IBM/flowdiscovery-simulator>. Code for running simulations using the Docker backend of ST4SD is available at <https://github.com/st4sd/flow-simulator-experiment>.

Related Publications
------------
* “Optimizing carbon dioxide trapping for geological storage.” Jaione Tirapu Azpiroz, et al., 2023. arXiv preprint arXiv:2312.13512, <https://doi.org/10.48550/arXiv.2312.13512>
* “Geometry evolution of porous media due to coupled reactive-transport processes within capillary networks.” David Alejandro Lazo Vasquez, et al. ACS Fall 2023. <https://doi.org/10.1021/scimeetings.3c10238>
* “Modeling carbon dioxide trapping at microscopic pore scale with digital rock representations.” Jaione Tirapu-Azpiroz et al. Proceedings of the SPIE, 12374-14, 2023. <https://doi.org/10.1117/12.2650243>
* “Full scale, microscopically resolved tomographies of sandstone and carbonate rocks augmented by experimental porosity and permeability values.” Esteves Ferreira, M., et al. Sci Data 10, 368 (2023). <https://doi.org/10.1038/s41597-023-02259-z>
* “Full scale, microscopically resolved tomographies of sandstone and carbonate rocks augmented by experimental porosity and permeability values.” Esteves Ferreira, Matheus, et al. (2022): Figshare+. Dataset. <https://doi.org/10.25452/figshare.plus.21375565.v5>
* “Cloud-based pore-scale simulator for studying carbon dioxide flow in digital rocks.” Tirapu Azpiroz, Jaione, et al. Proceedings of the 16th Greenhouse Gas Control Technologies Conference (GHGT-16) 23-24 Oct 2022,  <http://dx.doi.org/10.2139/ssrn.4276744>
* “Micro-computed tomography of sandstone rocks: Raw, filtered and segmented datasets” Everton Lucas-Oliveira, et al. Data in Brief, Volume 41, 2022, 107893, <https://doi.org/10.1016/j.dib.2022.107893>.
* "Advanced optical on-chip analysis of fluid flow for applications in carbon dioxide trapping," Jaione Tirapu-Azpiroz, et al., Proc. SPIE 11955, 1195507 (2022); <https://doi.org/10.1117/12.2610336>
* "High accuracy capillary network representation in digital rock reveals permeability scaling functions.” Neumann, R.F., et al. Sci Rep 11, 11370 (2021). <https://doi.org/10.1038/s41598-021-90090-0>
* "11 Sandstones: raw, filtered and segmented data", R. Neumann, et al. Digital Rocks Portal, 2020. <https://dx.doi.org/10.17612/f4h1-w124>
