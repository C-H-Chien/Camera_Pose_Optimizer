## Camera Pose Optimizer
### Introduction
This repository hosts a few camera pose optimizer for visual odometry or SLAM applications. They are isolated CPP functions with specified inputs and outputs for an easy-to-plug-in purpose, _e.g._, used as a wrapper function communicated with the python or the MATLAB code. 

### Dependencies
Eigen 3.X

### Optimization Solvers
**Semi-Direct Photometry-Based Solver**: See `optimize_pose_mex.cpp`. 

### Build up a mex file for use in MATLAB
For each solver to be used by MATLAB, build it in the MATLAB console by, for example, the semi-direct photometry-based solver,
```bash
mex solver_for_matlab/optimize_pose_mex.cpp solver_for_matlab/photometry_minimizer.cpp -I/path/of/eigen3/install/prefix
```
where `/path/of/eigen3/install/prefix` is the install prefix of the eigen library, _e.g._, `-I/usr/include/eigen3/`.

