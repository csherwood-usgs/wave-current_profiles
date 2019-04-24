# wave-current_profiles

### Near-bottom profiles generated from simple wave-current models
---
This resides in my ../src folder.

`m94_orbital.m` - This generates wave-orbital profiles for various parts of the wave cycle and mean current profile in the along-wave and cross-wave direction. The model (`Madsen94.pdf`) is based on a discountinuous eddy-viscosity profile, so there is an unrealistic amount of shear across the wave boundary layer. The code requires these functions: `m94r.m`, `ker.m`, and `kei.m`.

`smith77.m` - This is a script that produces profiles of wave-orbital velocities only based on Smith, J.D. (1977) Chap. 13 in The Sea, vol. 6. The complete model is not here, because you need to know the wave-current _u*_. Maybe TODO: add that calculation. 

`klvn.m' - Every time I mess with this, I worry about the Kelvin functions, so I keep this code as a touchstone.

TODO: add the Smith-Kachel-Wiberg model, which has continuous eddy viscosity and velocity profiles
