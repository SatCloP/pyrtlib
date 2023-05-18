<img align="" src="resources/logo/logo_large_new.png" width="400">

# A Radiative Transfer Python Library (non-scattering)

[![docker-image-ci](https://github.com/slarosa/pyrtlib/workflows/docker-image-ci/badge.svg)](https://github.com/slarosa/pyrtlib/actions/workflows/docker-image.yml)
[![run-python-tests](https://github.com/slarosa/pyrtlib/workflows/run-python-tests/badge.svg)](https://github.com/slarosa/pyrtlib/actions/workflows/ci.yml)
[![build-docs-action](https://github.com/slarosa/pyrtlib/workflows/build-docs-action/badge.svg)](https://github.com/slarosa/pyrtlib/actions/workflows/build_docs.yml)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/slarosa/pyrtlib?display_name=tag)
[![codecov](https://codecov.io/gh/slarosa/pyrtlib/branch/main/graph/badge.svg?token=7DV4B4U1OZ)](https://codecov.io/gh/slarosa/pyrtlib)
[![license](https://img.shields.io/github/license/slarosa/pyrtlib.svg)](https://github.com/slarosa/pyrtlib/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234.svg)](https://doi.org/10.5281/zenodo.1234)
<!--[![GitHub commit](https://img.shields.io/github/last-commit/slarosa/pyrtlib)](https://github.com/slarosa/pyrtlib/commits/main)-->
<!-- [![license](https://img.shields.io/github/license/slarosa/pyrtlib.svg)](https://github.com/slarosa/pyrtlib/blob/main/LICENSE.md) -->

A python package to compute atmospheric radiative transfer model based on the radiative transfer equation (RTE).

![spectrum](resources/spectrum_r22.jpeg)

# Example

For examples of how to use pyrtlib see the [examples gallery](docs/examples). Code can be downloaded both as python script or notebook file.

## Performing calculation of upwelling brightness temperature.

```python
   from pyrtlib.tb_spectrum import TbCloudRTE
   from pyrtlib.atmospheric_profiles import AtmosphericProfiles as atmp
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
