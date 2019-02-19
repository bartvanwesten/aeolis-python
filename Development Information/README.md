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
- Vegetation establishment (germination + distribution)
- Vegetation growth (Duran, 2006)
- Vegetation lateral propagation (DUBEVEG)
- Vegetation shear stress reduction (Raupach, 1993)
- Influence of hydrodynamics on morphology and vegetation

These implementations were mainly derived from CDM (Coastal Dune Model) and DUBEVEG (DUne-BEach-VEGetation).

The new AeoLiS version is now capable of simulating several aspects of dune development, 
including the formation of barchan dunes, parabolic and embryo dune fields.

# Installation

### Prerequisites 
1. Install Anaconda from https://www.continuum.io/downloads (Python 3.x version).
2. Create an environment 'aeolis' in the Anaconda Navigator.
3. Open the Anaconda Prompt and run 'activate aeolis' to activate the enviroment.
4. Run the following commands:

```
 pip install numpy
 pip install scipy
 pip install netcdf4
 pip install bmi-python
```
  
### Installation of AeoLiS
5. The installation-method of AeoLiS depends on the desired AeoLiS-version

5a. "Original" (Hoonhout, 2017) AeoLiS-version. In Anaconda Prompt run:
```
pip install aeolis
```

5b. Alternative AeoLiS (newer) version. Download/clone the desired AeoLiS version and in Anaconda Prompt run:

Development
```
python setup.py develop
```
Static installation
```
python setup.py install
```
    
# Model Limitations and possible solutions
Since the model is still in development (and will be for a long time) it has its limitations.
These limitations can be divided into two categories" (for more information, see chapter 5-7 in http://resolver.tudelft.nl/uuid:83344d5b-f232-4085-985c-d1ec38903254):

### "Academic/simplistic" dune development
The current implementation of processes had the aim to make the model capable of simulating academic dune development (barchan, parabolic, embryo dunes). In order to validate the simulation results, results from CDM are used. The initial simulation results are reasonably well, but differ slightly from CDM and computations are significantly slower. In order to improve this, the following aspects of the model should be considered:
  - Advection equation
  - Multi-fraction sediment transport
  - Grid resolution
  - Computational speed (10x-100x, not only workability also increases work-range)
  
### "Realistic/complex" dune development
AeoLiS has the potential to simulate the behaviour of complete coastal areas, but it still lacks the description of many processes. In order to increase the quantitative predicitability of the model, the implementation of the following processes should be consired. More information on these processes can be found in chapter 7 of http://resolver.tudelft.nl/uuid:83344d5b-f232-4085-985c-d1ec38903254. 

# List of relevant literature

### CDM

Durán, O. (2007). Vegetated dunes and barchan dune fields.

Parteli, E. J., Kroy, K., Tsoar, H., Andrade, J. S., & Pöschel, T. (2014). Morphodynamic modeling of aeolian dunes: Review and future plans. The European Physical Journal Special Topics, 223(11), 2269-2283.

### AeoLiS

Hoonhout, B. (2017). Aeolian Sediment Availability and Transport. https://doi.org/10.4233/uuid:e84894d6-87d2-4006-a8c2-d9fbfacabddc

van Westen, B. (2018). “Numerical modelling of aeolian coastal landform development,” M.Sc. Thesis. Technical University of Delft. Delft, the Netherlands. http://resolver.tudelft.nl/uuid:83344d5b-f232-4085-985c-d1ec38903254 

### DUBEVEG

De Groot, A. V., Berendse, F., Riksen, M., Baas, A., Slim, P. A., Van Dobben, H. F., & Stroosnijder, L. (2011). Modelling coastal dune formation and associated vegetation development.

### Background of physical processes and implementations

Durán, O., & Herrmann, H. J. (2006). Vegetation against dune mobility. Physical review letters, 97(18), 188001.

Raupach, M. R., Gillette, D. A., & Leys, J. F. (1993). The effect of roughness elements on wind erosion threshold. Journal of Geophysical Research: Atmospheres, 98(D2), 3023-3029.

Sauermann, G., Kroy, K., & Herrmann, H. J. (2001). Continuum saltation model for sand dunes. Physical Review E, 64(3), 031305.

Weng, W. S., Hunt, J. C. R., Carruthers, D. J., Warren, A., Wiggs, G. F. S., Livingstone, I., & Castro, I. (1991). Air flow and sand transport over sand-dunes. In Aeolian Grain Transport (pp. 1-22). Springer, Vienna.


