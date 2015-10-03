# Octopus scripts

A collection of scripts to postprocess data generated whith Octopus (www.octopus-code.org), a highly scalable rea-space real-time TDDFT code.


## bands2h.py
This utility converts the bands structure file (```bands-gp.dat```) generated by Octopus into a human-readable and easy-to-plot format. The script is capable to handle 1D, 2D and 3D data and perform cuts interpolating the band-structure on arbitrary segmented-lines in the Brillouin zone. 

#### Examples
To cut the band structure contained in ```file``` along the segment connecting (0,0,0) and (0,0,1/2) in reduced coordinates:  
```
bands2h.py -r -c "0,0,0; 0,0,1/2" file
```  

To generate the complete band structure from 2D data:
```
bands2h.py file
```  
 

For more info and a complete list of features:  
```
bands2h.py --help
``` 

#### External dependencies
* [numpy](http://www.numpy.org)
* [scipy](http://www.scipy.org)

## wfcutter.py
This utility allow to analyze wavefunctions and density obtained with octopus integrating the charge over volumes defined as the union of basic shapes. It also provides a basic interface to visualize the data and the integration volume. 

NOTE: It only supports vtk data format.

#### Examples
To integrate a real wavefunction contained in ```file.vtk``` over a sphere of radius 10 centered in (0,0,0):
```
wfcutter.py -g "sphere(0,0,0,10)" file.vtk

```

To visualize an isosurface of the data and the integration volume (defined by two interloking spheres):
```
wfcutter.py wfcutter.py -g "sph(0,0,0,10); sph(5,0,0,10)" -d file.vtk

```

For more info and a complete list of features:  
```
wfcutter.py --help
``` 

#### External dependencies
* [numpy](http://www.numpy.org)
* [scipy](http://www.scipy.org)
* [vtk](http://www.vtk.org/Wiki/VTK)
