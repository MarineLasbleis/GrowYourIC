# GrowYourIC
From geodynamics to Seismic observations in the Earth's inner core

![RotationTranslationGrowth](https://github.com/MarineLasbleis/GrowYourIC/blob/master/RTP.png "RotationTranslationGrowth")


This program calculates seismic observations (synthetics one) from
a given geodynamical model.

Calculate raypath:
- from real data set (data file has to be given)
 * from Lauren's data set : station, PKIKP-PKiKP travel time residual, zeta (angle from rotation axis), epicentral distance, station lat, station lon, event lat, event lon, event depth, inner core in lat, in lon, out lat, out lon, turn lat, turn lon, turn depth (in inner core), inner core travel time, PKIKP/PKiKP amplitude ratio
- for random sampling (surface and depth repartition can be choosed or
set to default values)
- for perfect sampling (either only in the equatorial plane, or in the total volume)
- TO DO: perfect sampling but with Poisson-disc sampling (see Mitchell’s Best-Candidate)

Geodynamical models:
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



files:
- positions.py
 * Point()
  + SeismoPoint(Point)
  + CartesianPoint(Point)
  + RandomPoint(Point)
	Raypath()
		Raypath_BT(Raypath)
		Raypath_inout(Raypath)
- geodynamic.py
	ModelGeodynamic()
		PureTranslation(ModelGeodynamic)
		TranslationRotation(ModelGeodynamic)
		PureRotation(ModelGeodynamic)
		PureGrowth(ModelGeodynamic)
		TranslationGrowth(ModelGeodynamic)
		TranslationGrowthRotation(ModelGeodynamic)
- data.py
	SeismicData()
		SeismicFromFile(SeismicData)
		PerfectSamplingEquator(SeismicData)
		RandomData(SeismicData)
- intersection.py



 TO DO:

- add other proxies and way to compute proxies. Now, only age is rendered. (geodynamic.py)
- add other intersections methods. Only brentq is used. see intersection.py. Also, adding a way to remember the choice of the user could be a good idea (for the values of the limits of the brentq) For example, first plot the trajectory_r compared to radius for N (10? 100?) samples, and allow the user to input the [x0, x1] he want to use for the run. 
