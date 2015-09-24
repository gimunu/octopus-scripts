# Octopus scripts

A collection of scripts to postprocess data generated from Octopus (www.octopus-code.org), an highly scalable rea-space real-time TDDFT code.


## bands2h.py
This utility converts the bands structure file generated by Octopus into a human-readable and easy-to-plot format. The script is capable to handle 1D, 2D and 3D data and perform cuts interpolating the band-structure on arbitrary segmented-lines in the Brillouin zone. 

#### Dependencies
The script needs [numpy](http://www.numpy.org) and [scipy](http://www.scipy.org) to run.
