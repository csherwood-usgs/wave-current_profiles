# wave-current_profiles

### Near-bottom profiles generated from simple wave-current models
This resides in my ../src folder.

---

`m94_orbital.m` - This generates wave-orbital profiles for various parts of the wave cycle and mean current profile in the along-wave and cross-wave direction. The model (`Madsen94.pdf`) is based on a discountinuous eddy-viscosity profile, so there is an unrealistic amount of shear across the wave boundary layer. The code requires these functions: `m94r.m`, `ker.m`, and `kei.m`. In this example, current velocity measured 1 m above the bed is 0.2 m/s, and the wave-orbital velocity is 0.5 m/s. The various colored lines are orbital-velocity profiles at various points in the wave cycle, the dark blue line is the mean current profile in the direction of the waves, and the red line (zeros) is the mean current perpendicular to the wave direction. The waves and currents are aligned, so the red profile is all zeros. The top of the wave boundary layer is shown with a dashed line.

![alt text](https://github.com/csherwood-usgs/wave-current_profiles/raw/master/profs.png "Profiles")

Here are the profiles (along- and across-wave directions) at point in the wave cycle where waves and currents combine to create the maximum speeds.

![alt text](https://github.com/csherwood-usgs/wave-current_profiles/raw/master/comb_profs.png "Max. Profiles")


---

`smith77.m` - This is a script that produces profiles of wave-orbital velocities only based on Smith, J.D. (1977) Chap. 13 in The Sea, vol. 6. The complete model is not here, because you need to know the wave-current _u*_. Maybe TODO: add that calculation. 

---

`klvn.m' - Every time I mess with this, I worry about the Kelvin functions, so I keep this code as a touchstone.

---

TODO: add the Smith-Kachel-Wiberg model, which has continuous eddy viscosity and velocity profiles

---
