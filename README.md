<img align="" src="https://raw.githubusercontent.com/SatCloP/pyrtlib/main/resources/logo/logo_large_new.png" width="400">

# A Radiative Transfer Python Library (non-scattering)

[![docker-image-ci](https://github.com/SatCloP/pyrtlib/actions/workflows/docker-image.yml/badge.svg)](https://github.com/SatCloP/pyrtlib/actions/workflows/docker-image.yml)
[![run-python-tests](https://github.com/SatCloP/pyrtlib/actions/workflows/ci.yml/badge.svg)](https://github.com/SatCloP/pyrtlib/actions/workflows/ci.yml)
[![build-docs-action](https://github.com/SatCloP/pyrtlib/actions/workflows/build_docs.yml/badge.svg)](https://github.com/SatCloP/pyrtlib/actions/workflows/build_docs.yml)

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![PyPI Latest Release](https://img.shields.io/pypi/v/pyrtlib.svg)](https://pypi.org/project/pyrtlib/)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/SatCloP/pyrtlib?display_name=tag)

[![codecov](https://codecov.io/gh/SatCloP/pyrtlib/branch/main/graph/badge.svg?token=7DV4B4U1OZ)](https://codecov.io/gh/SatCloP/pyrtlib)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/55e8f5d14a83477aaab2eff090b10281)](https://app.codacy.com/gh/SatCloP/pyrtlib/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

[![license](https://img.shields.io/github/license/SatCloP/pyrtlib.svg)](https://github.com/SatCloP/pyrtlib/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/345925671.svg)](https://zenodo.org/badge/latestdoi/345925671)

<!--[![GitHub commits since tagged version](https://img.shields.io/github/commits-since/SatCloP/pyrtlib/v1.0.0)](https://github.com/SatCloP/pyrtlib/commits/) -->
<!--[![GitHub commit](https://img.shields.io/github/last-commit/slarosa/pyrtlib)](https://github.com/SatCloP/pyrtlib/commits/main)-->
<!-- [![license](https://img.shields.io/github/license/slarosa/pyrtlib.svg)](https://github.com/SatCloP/pyrtlib/blob/main/LICENSE.md) -->

PyRTlib is a Python package, for non-scattering line-by-line microwave RT simulations. PyRTlib is a user-friendly tool for computing down and up-welling brightness temperatures and related quantities (e.g., atmospheric absorption, optical depth, opacity) in Python.

![spectrum](https://raw.githubusercontent.com/SatCloP/pyrtlib/main/resources/spectrum_r22.jpeg)

Plotting of nadir upwelling $\Delta T_b$ using the last two absorption models available in PyRTlib for six reference atmosphere climatology.

![spectrum](https://raw.githubusercontent.com/SatCloP/pyrtlib/main/resources/spectrum_r23_r24.png)


# Installation

Use pip package to install quicly the pyrtlib library. See [installation instructions](https://satclop.github.io/pyrtlib/en/main/installation.html) for more info on how instaling pyrtlib. 

```sh
   $ pip install pyrtlib
```

# Example

For examples of how to use pyrtlib see the [examples gallery](https://satclop.github.io/pyrtlib/en/main/examples/index.html). Code can be downloaded both as python script or notebook file.

## Performing calculation of upwelling brightness temperature.

```python
   from pyrtlib.tb_spectrum import TbCloudRTE
   from pyrtlib.climatology import AtmosphericProfiles as atmp
   from pyrtlib.utils import mr2rh, ppmv2gkg
```
Atmospheric profile definition:

```python   
   z, p, _, t, md = atmp.gl_atm(atmp.MIDLATITUDE_SUMMER)
```

Units conversion:

```python 
   gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
```
Relative humidity of $H_2O$ (water vapor)

```python
   rh = mr2rh(p, t, gkg)[0] / 100
```
Deifinition of angles and frequencies:

```python 
   ang = np.array([90.])
   frq = np.arange(20, 1001, 1)
```
Initialize parameters for main execution:

```python 
   rte = TbCloudRTE(z, p, t, rh, frq, ang)
```
Set absorption model:

```python 
   rte.init_absmdl('R22SD')
```
Execute model by computing upwelling radiances:

```python 
   df = rte.execute()
   df.tbtotal
   0      293.119811
   1      292.538088
   2      291.736672
   3      291.913658
   4      292.493971
            ...    
   976    230.179993
   977    231.435965
   978    232.592915
   979    233.666322
   980    234.667522
   Name: tbtotal, Length: 981, dtype: float64
```

## My first run with PyRTlib

You can get started with PyRTlib by installing and executing the first radiative transfer calculation from the following Colab Notebook [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SatCloP/pyrtlib/blob/main/docs/source/notebook/first_run.ipynb)


## Cite as

Larosa, S., Cimini, D., Gallucci, D., Nilo, S. T., and Romano, F.: PyRTlib: an educational Python-based library for non-scattering atmospheric microwave radiative transfer computations, Geosci. Model Dev., 17, 2053–2076, https://doi.org/10.5194/gmd-17-2053-2024, 2024.

Larosa, S., Cimini, D., Gallucci, D., Nilo, S. T., & Romano, F. (2024). PyRTlib: a python package for non-scattering line-by-line microwave Radiative Transfer simulations. (Computer software). https://doi.org/10.5281/zenodo.8219145

## Contributors
<a href="https://github.com/SatCloP/pyrtlib/graphs/contributors"><img align="" src="https://contrib.rocks/image?repo=SatCloP/pyrtlib"></a>
