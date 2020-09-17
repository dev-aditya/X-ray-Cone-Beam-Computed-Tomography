# X-ray Cone Beam Computed Tomography

X-ray Cone-beam Computed tomography (CBCT) is a technique for imaging cross-sections of an object using a series of X-ray measurements taken from
different angles around the object. CBCT uses a cone shaped X-ray beam.
The projection data acquired at the 2D detector array, which is then reconstructed into a 3D image volume using the FDK reconstruction
algorithm.

This repository contains computer program to simulate CBCT. 


```
├── README.md
├── LICENSE
├── MATLAB
│   ├── main
│   │   ├── back_project.m
│   │   ├── phantom.m
│   │   ├── projections.m
│   │   ├── ramp_filter.m
│   │   ├── shepp_logan.m
│   │   └── ye_yu_wang.m
│   ├── parameter.mat
│   ├── xray_cone_beam_reconstruction_filtered.mlx
│   └── xray_cone_beam_reconstruction.mlx
└── Python
    ├── back_projection.py
    ├── phantom_features.py
    ├── phatom_const.py
    ├── projection.py
    └── ramp.py
```
