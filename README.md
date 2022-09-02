# HAYASHI

### (Halo-level AnalYsis of the Absorption Signal in HI)

(hayashi means *forest* in japanese)

[![arXiv](https://img.shields.io/badge/arXiv-22XX.XXXXX-B31B1B.svg)](http://arxiv.org/abs/22XX.XXXXX) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4707447.svg)](https://doi.org/10.5281/zenodo.4707447)

Python library for computing the number of absorption features of the 21 cm forest in a semianalytic formalism. Includes the enhancement of the signal due to the presence of substructures within minihalos, as studied in [arXiv:22XX.XXXXX](https://arxiv.org/abs/22XX.XXXXX), see that paper for more details.

## Notebook examples

In order to illustrate the usage of the library, we include several example notebooks:

* `absorbers_example.ipynb`: computes all the relevant outputs for the 21 cm forest, such as the optical depth, the maximum impact parameter and the number of absorbers, comparing the cases with and without the subhalo contribution.

* `tidal_disruption.ipynb`: compares the 21 cm forest observables when tidal disruption is considered in subhalos.

* `density_profiles.ipynb`: compares the 21 cm forest outputs for different density profiles: NFW and uniform.

## Requisites

The code is written in Python3, and makes use of the package for cosmological computations [Colossus](https://bdiemer.bitbucket.io/colossus/), as well as several Python libraries:

* `numpy`
* `matplotlib`
* `scipy`
* `tqdm`


## Citation

If you use the code, please link this repository, and cite [arXiv:22XX.XXXXX](https://arxiv.org/abs/22XX.XXXXX) and the DOI [10.5281/zenodo.4707447](https://doi.org/10.5281/zenodo.4707447).

## Contact

For comments, questions etc. you can contact me at <pablo.villanueva.domingo@gmail.com>
