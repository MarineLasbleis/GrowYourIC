

[![DOI](https://zenodo.org/badge/49926307.svg)](https://zenodo.org/badge/latestdoi/49926307)


# GrowYourIC
From geodynamics to Seismic observations in the Earth's inner core

![RotationTranslationGrowth](https://github.com/MarineLasbleis/GrowYourIC/blob/master/RTP.png "RotationTranslationGrowth")


Associated publication (please cite it if you are using this toolbox for scientific purpose)
M. Lasbleis, L. Waszek, E. Day. GrowYourIC: a step towards a coherent model of seismic structure,
accepted in Geochemistry, Geophysics, Geosystems. DOI:10.1002/2017GC007149


## Overview
This program calculates seismic observations (synthetics one) from
a given geodynamical model.

Calculate raypath (only PKIKP types paths are implemented now):
- from real data set (data file has to be given)
 * from Lauren's data set : station, PKIKP-PKiKP travel time residual, zeta (angle from rotation axis), epicentral distance, station lat, station lon, event lat, event lon, event depth, inner core in lat, in lon, out lat, out lon, turn lat, turn lon, turn depth (in inner core), inner core travel time, PKIKP/PKiKP amplitude ratio
- for random sampling (surface and depth repartition can be choosed or
set to default values)
- for perfect sampling (either only in the equatorial plane, or in the total volume). "perfect" sampling is used to construct meshgrids and contourf plots. 
- TO DO: perfect sampling but with Poisson-disc sampling (see Mitchell’s Best-Candidate)

Geodynamical models:
- Translation (Monnereau et al. 2010, Alboussiere at al. 2010, Geballe
et al. 2013, Deguen et al. 2014) (return age or growth rate of material)
- Translation - growth - super rotation (return age or growth rate of material)

Mineral Physics:
- Calvet and Margerin 2008 data (polynomial fit of the figures of the paper)


For now, the code does not provide "seismic observations" (=seismic velocity and/or seismic residuals), 
but we try
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

## Quickstart

### Python dependencies:
(please use pip to install them)

numpy

scipy

matplotlib (version 2.5 at least)

matplotlib.basemap (to be installed from their website : http://matplotlib.org/basemap/ 

### installation: 
>> python3 setup.py install

to check the installation:
python3 examples.py
(does not require external database for seismic rays)

### Data distributed with the package

The folder GrowYourIC/data contains a couple of data files that can be used by the user 
(if used for any type of publication, please refer to each file header with citation information. 
Do not distribute the files without their header and proper citation.)

Please copy the data files in your work folder to be able to use them. Examples may be run in the folder GrowYourIC and will use directly data files in the original folder. 

## List of files

files: (classes start with a capital letter, regular functions with a small letter)
- positions.py
 + from_seismo_to_cartesian()
 + from_cartesian_to_seismo()
 + angular_distance_to_point()
 + Point(object)
 + SeismoPoint(Point)
 + CartesianPoint(Point)
 + RandomPoint(Point)
 + Raypath(object)
 + Raypath_BT(Raypath)
 + Raypath_inout(Raypath)
- geodyn.py
 + evaluate_proxy()
 + average_proxy()
 + Model(object)
- geodyn_trg.py
 + translation_velocity()
 + radial_derivative()
 + ModelTRG(geodyn.Model)
 + PureTranslation(ModelTRG)
 + TranslationRotation(ModelTRG)
 + PureRotation(ModelTRG)
 + PureGrowth(ModelTRG)
 + TranslationGrowth(ModelTRG)
 + TranslationGrowthRotation(ModelTRG)
- geodyn_static.py
 + Hemispheres(geodyn.Model)
- data.py
 + read_from_file()
 + SeismicData(object)
 + SeismicFromFile(SeismicData)
 + PerfectSamplingEquator(SeismicData)
 + RandomData(SeismicData)
 + PerfectSamplingEquatorRadial(SeismicData)
- intersection.py
 + zero_brentq() (only function used outside the module)
- mineral_phys_data.py
 + export_matlab_data(name_data, file_name="./CM2008/data.mat"):
 + domain_size(age):
 + adimensional_frequency(size, v=11030., freq=1.):
 + convert_CM2008_velocity(kR, poly):
 + convert_CM2008_attenuation(kR, poly):
 + heaviside(x):
- plot_data.py
 + setting_map() #to prepare a map for plotting data.


 TO DO:

- add other intersections methods. Only brentq is used. see intersection.py. Also, adding a way to remember the choice of the user could be a good idea (for the values of the limits of the brentq) For example, first plot the trajectory_r compared to radius for N (10? 100?) samples, and allow the user to input the [x0, x1] he want to use for the run. 
