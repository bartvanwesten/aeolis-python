[![CircleCI](https://circleci.com/gh/openearth/aeolis-python.svg?style=svg)](https://circleci.com/gh/openearth/aeolis-python)
[![Codecov](https://codecov.io/gh/openearth/aeolis-python/branch/master/graph/badge.svg)](https://codecov.io/gh/openearth/aeolis-python)
[![ReadTheDocs](http://readthedocs.org/projects/aeolis/badge/?version=latest)](http://aeolis.readthedocs.io/en/latest/)

[![PyPI](https://img.shields.io/pypi/v/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![PyPI_versions](https://img.shields.io/pypi/pyversions/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![PyPI_status](https://img.shields.io/pypi/status/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![PyPI_format](https://img.shields.io/pypi/format/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)

[![License](https://img.shields.io/pypi/l/AeoLiS.svg)](https://pypi.python.org/pypi/AeoLiS)
[![DOI](https://zenodo.org/badge/7830/openearth/aeolis-python.svg)](https://zenodo.org/badge/latestdoi/7830/openearth/aeolis-python)

# AeoLiS
AeoLiS is a process-based model for simulating aeolian sediment
transport in situations where supply-limiting factors are important,
like in coastal environments. Supply-limitations currently supported
are soil moisture contents, sediment sorting and armouring and
roughness elements.

Documentation can be found at
http://aeolis.readthedocs.io/.

AeoLiS is initially developed by [Bas Hoonhout](mailto:b.m.hoonhout@tudelft.nl)
at Delft University of Technology with support from the ERC-Advanced
Grant 291206 Nearshore Monitoring and Modeling
([NEMO](http://nemo.citg.tudelft.nl>)) and
[Deltares](http://www.deltares.nl>). AeoLiS is currently maintained by
[Bas Hoonhout](mailto:bas.hoonhout@deltares.nl) at Deltares and
[Sierd de Vries](mailto:Sierd.deVries@tudelft.nl) at Delft University of Technology.

## Examples

```
aeolis params.txt
```

# AeoLiS development

The following physical processes have been recently implemented:
- Saltation model (Sauermann, 2001)
- Analytical perturbation theory (Weng, 1991)
- Separation bubble (Sauermann, 2001)
- User-defined non-erodible layer
- Avalanching
- 

These implementations mainly were derived from CDM (Coastal Dune Model) and DUBEVEG (DUne-BEach-VEGetation).

The new AeoLiS version is now capable of simulating several aspects of dune development, 
including the development of barchan dunes, parabolic and embryo dune fields.

# Installation

## Prerequisites 

1. Install Anaconda from https://www.continuum.io/downloads (Python 3.x version).
2. Create an environment 'aeolis' in the Anaconda Navigator.
3. Open the Anaconda Prompt and run 'activate aeolis' to activate the enviroment.
4. Run the following commands:
  - pip install numpy
  - pip install scipy
  - pip install netcdf4
  - pip install bmi-python
  
5. The installation-method of AeoLiS depends on the desired AeoLiS-version
