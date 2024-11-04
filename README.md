# HAYASHI

### (Halo-level AnalYsis of the Absorption Signal in HI)

[![arXiv](https://img.shields.io/badge/arXiv-2209.01305-B31B1B.svg)](http://arxiv.org/abs/2209.01305) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7044255.svg)](https://doi.org/10.5281/zenodo.7044255) ![Python](https://img.shields.io/pypi/pyversions/python-binance.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

(hayashi means *forest* in japanese)

Python library for computing the number of absorption features of the 21 cm forest in a semianalytic formalism. Includes the enhancement of the signal due to the presence of substructures within minihalos, as studied in [arXiv:2209.01305](https://arxiv.org/abs/2209.01305). It supports non-standard cosmologies with impact in the large scale structure, such as warm dark matter and primordial black holes. See the papers [arXiv:2209.01305](https://arxiv.org/abs/2209.01305), [arXiv:2104.10695](https://arxiv.org/abs/2104.10695) for more details.

[Read the documentation here.](https://hayashi.readthedocs.io/en/latest/)


## Installation

The code is written in Python3, and makes use of the package for cosmological computations [Colossus](https://bdiemer.bitbucket.io/colossus/), as well as several standard Python libraries(`numpy`,`scipy`,`tqdm`), which are automatically installed when `hayashi` is installed.

For installing the Python package from [PyPI](https://pypi.org/project/hayashi/):

`pip install hayashi`


## Usage

The basis of the code is the 21 cm `Forest` class. Given a redshift and the temperature of the intergalactic medium at that epoch, we can define an instance of the state of the 21 cm forest.

```python
from hayashi.forest import Forest
from hayashi.cosmo import Tk_ad

# Define the redshift of interest
z = 10
# Get the adiabatic temperature of the intergalactic medium at z
Tk = Tk_ad(z)

21cmforest = Forest(z, Tk)
```

This allows to call different observables such as the optical depth or the number of absorbers.

```python
# Get a the optical depth, as a matrix in mass and impact parameter
tau = 21cmforest.tau_tot

# Get the number of absorption features and its (logarithmic) derivative with respect to tau
Nabs, dNabsdtau = 21cmforest.num_absorbers()
```

It is straightforward to include non-standard cosmologies by replacing the halo mass function, either using those included in the code or defined by the user. This is an example with primordial black holes, which modify the halo mass function due to a shot noise isocurvature mode (see [arXiv:2104.10695](https://arxiv.org/abs/2104.10695)):

```python
from hayashi.nlcdm import dndlnM_PBH

# Define a cosmology where 10 % of dark matter is composed by primordial black holes of 1 solar mass
21cmforest_PBH = Forest(z, Tk, dndlnM = lambda M, z: dndlnM_PBH(M, z, fpbh = 0.1, Mpbh = 1.))
```

See the source code at `hayashi` for more details, and the sample notebooks for examples of usage.


## Notebook examples

In order to illustrate the usage of the library, we include several example notebooks:

* `absorbers_example.ipynb`: computes all the relevant outputs for the 21 cm forest, such as the optical depth, the maximum impact parameter and the number of absorbers, comparing the cases with and without the subhalo contribution.

* `nlcdm_example.ipynb`: compares the standard CDM 21 cm forest with different non-standard cosmologies: warm dark matter and primordial black holes.

* `density_profiles.ipynb`: compares the 21 cm forest outputs for different density profiles: NFW and uniform.

* `tidal_disruption.ipynb`: compares the 21 cm forest observables when tidal disruption is considered in subhalos.


## Citation

If you use the code, please link this repository, and cite [arXiv:2209.01305](https://arxiv.org/abs/2209.01305) and the DOI [10.5281/zenodo.7044255](https://doi.org/10.5281/zenodo.7044255).


## Contact

For comments, questions etc. you can contact me at <pablo.villanueva.domingo@gmail.com>
