# Non-ΛCDM models

It is straightforward to include non-standard cosmologies by replacing the halo mass function, either using those included in the code or defined by the user.

The current implementation includes two different non-ΛCDM scenarios:
- Warm Dark Matter (WDM)
- Primordial Black Holes (PBH)

This is an example with primordial black holes, which modify the halo mass function due to a shot noise isocurvature mode (see [arXiv:2104.10695](https://arxiv.org/abs/2104.10695)):

```py
from hayashi.nlcdm import dndlnM_PBH

# Define a cosmology where 10 % of dark matter is composed by primordial black holes of 1 solar mass
21cmforest_PBH = Forest(z, Tk, dndlnM = lambda M, z: dndlnM_PBH(M, z, fpbh = 0.1, Mpbh = 1.))
```

See the source code at `hayashi` for more details, and the sample notebooks for examples of usage.
