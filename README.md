<img align="" src="resources/logo/logo_large_new.png" width="400">

# A Radiative Transfer Python Library (non-scattering)

[![build-docs-action](https://github.com/slarosa/pyrtlib/workflows/build-docs-action/badge.svg)](https://github.com/slarosa/pyrtlib/actions/workflows/build_docs.yml)
[![run-python-tests](https://github.com/slarosa/pyrtlib/workflows/run-python-tests/badge.svg)](https://github.com/slarosa/pyrtlib/actions/workflows/ci.yml)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![GitHub release](https://img.shields.io/github/release/slarosa/pyrtlib.svg)](https://github.com/slarosa/pyrtlib)
[![codecov](https://codecov.io/gh/slarosa/pyrtlib/branch/main/graph/badge.svg?token=7DV4B4U1OZ)](https://codecov.io/gh/slarosa/pyrtlib)
<!--[![GitHub commit](https://img.shields.io/github/last-commit/slarosa/pyrtlib)](https://github.com/slarosa/pyrtlib/commits/main)-->
<!-- [![license](https://img.shields.io/github/license/slarosa/pyrtlib.svg)](https://github.com/slarosa/pyrtlib/blob/main/LICENSE.md) -->

A python package to compute atmospheric radiative transfer model based on the radiative transfer equation (RTE).

![spectrum](resources/spectrum.png)

# Example

For examples of how to use pyrtlib see the [examples gallery](docs/examples). Code can be downloaded both as python script or notebook file.

## Performing calculation of brightness temperature for a model.
Atmospheric profile definition:

```python
>>> z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)
```
Units conversion:
```python
>>> gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
```
Relative humidity of H2O (water vapor)
```python
>>> rh = mr2rh(p, t, gkg)[0] / 100
```
Deifinition of angles and frequencies:
```python
>>> ang = np.array([90.])
>>> frq = np.arange(20, 201, 1)
```
Initialize parameters:
```python
>>> rte = TbCloudRTE(z, p, t, rh, frq, ang)
```
Set absorption mofel:
```python
>>> rte.init_absmdl('R16')
```
Execute model:
```python
>>> df = rte.execute()
>>> df.tbtotal
0      297.391838
1      296.186240
2      294.748245
3      294.953483
4      296.027799
         ...
176    275.997899
177    276.611319
178    277.129218
179    277.566840
180    277.936645
Name: tbtotal, Length: 181, dtype: float64
```
Preview of the output dataframe:
|tbtotal|tbatm|tmr|tmrcld|tauwet|taudry|tauliq|tauice|
|:----|:----|:----|:----|:----|:----|:----|:----|
|20|297.391838|0.0|281.191287|0.0|0.120341|0.012855|0.0|0.0|
|21|296.186240|0.0|280.517137|0.0|0.188802|0.013524|0.0|0.0|
|22|294.748245|0.0|279.175653|0.0|0.261841|0.014259|0.0|0.0|
|23|294.953483|0.0|279.830575|0.0|0.257906|0.015066|0.0|0.0|
|24|296.027799|0.0|280.971991|0.0|0.202303|0.015954|0.0|0.0|
|...|...|...|...|...|...|...|...|...|
|196|275.997899|0.0|275.396235|0.0|3.672911|0.025784|0.0|0.0|
|197|276.611319|0.0|275.881854|0.0|3.459942|0.025956|0.0|0.0|
|198|277.129218|0.0|276.279020|0.0|3.289797|0.026129|0.0|0.0|
|199|277.566840|0.0|276.605436|0.0|3.152663|0.026302|0.0|0.0|
|200|277.936645|0.0|276.874793|0.0|3.041382|0.026476|0.0|0.0|
