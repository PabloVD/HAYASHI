# Usage

The basis of the code is the 21 cm `Forest` class. Given a redshift and the temperature of the intergalactic medium at that epoch, we can define an instance of the state of the 21 cm forest.

```python
from Source.forest import Forest
from Source.cosmo import Tk_ad

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