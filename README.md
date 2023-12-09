# Calculation of transmission coefficient
• The data obtained from the measurement is stored in a .dat file, which can be converted to a .root file using the convert.C macro or shell script convert. If there are
multiple files to be converted, the convert shell script can be run in the directory
where the files are located. All .dat files will be converted to .root files with the
suffix converted.root. <br>
• The transmission coefficient can be determined using the transmittance.C macro
or transmittance shell script. An output .root file containing the mean value (transmission coefficient), mean value error, and wavelength of the laser will be created.<br>
• To obtain a graph of the coefficient of total internal reflection as a function of
wavelength, where the results from measurements are compared to predictions
from scalar scattering theory, the transmittance graph.C macro can be used.
If there are multiple measurements, such as results from six different laser intensities, the hadd command can be applied to merge all six root files or the transmittance shell script can be used.
