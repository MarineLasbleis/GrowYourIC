# GrowYourIC
From geodynamics to Seismic observations in the Earth's inner core

This program calculates seismic observations (synthetics one) from
a given geodynamical model.

Calculate raypath:
- from real data set (data file has to be given)
- for random sampling (surface and depth repartition can be choosed or
set to default values)

Geodynamical model:
- Translation (Monnereau et al. 2010, Alboussiere at al. 2010, Geballe
et al. 2013, Deguen et al. 2014) (return age or growth rate of material)
- Translation - growth - super rotation (return age or growth rate of material)

Mineral Physics:
not yet implemented


For now, the code does not provide "seismic observations", but we try
to match the seismic observations with the spatial repartition of a
proxy which has a physical meaning (can be: age of material,
deformation of material, etc.)
The simplest example would be:
with the translation of the inner core
(Monnereau et al. 2010, Alboussiere et al. 2010), the age of material
inside the inner core can be mapped inside the IC; we use this age as
a proxy for the actual velocity of P-waves inside the IC, and argues
that isotropic P-waves velocity should have the same pattern as the IC
age to validate the translation hypothesis. We compute IC age average
over raypaths calculated from real data set and compare the pattern
with actual residual travel time from same data set. Next step would
be to predict growth rate of crystals at IC conditions and compute
P-waves travel time residuals (as done in Geballe 2013)






