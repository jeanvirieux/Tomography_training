# Tomography training based on ray approach and eikonal approach

Ray approach
============

Toy 2D codes for solving ordinary differential equations of ray tracing
written in Python3.

Simple analytical cases:
=======================

comparison between analytical solution and numerical solution based
on Runge-Kutta 2nd order scheme

>>> case of the model where the square of slowness has an arbitrary
constant gradient.


>>> case of the model where the velocity has a constant vertical gradient.


Ray tracing with shooting conditions:
====================================

Illustration of a numerical ray tracing algorithm based on Runge-Kutta solver
and a grid description of the velocity. The velocity is interpolated by
bspline functions. Such an interpolation allows the computation of spatial
first-derivatives (and also second-derivatives)


Ray tracing with boundary conditions:
====================================

Illustration of a numerical ray and paraxial ray tracing algorithm based on
Runge-Kutta solver and a grid description of the velocity. Finding initial
conditions at the source in order to reach the receiver is performed through
a Newton method where paraxial quantities provides the angle gradient with
respect to the ray-receiver distance.


Eikonal approach
================

Many codes exist for solving the Eikonal equation

The Eikonal equation can be solved efficiently for first-arrival times. Later
arrivals could be considered, but efficient approaches have to be found.

Fast marching approach is popular because the underlying causality implies
that the model inside which the time is computed is sampled only once

Fast sweeping approach is popular because there is no need to keep track of
the causal front evolution at the expense of iterative model sampling.

An efficient code is available in C++ by Jean-Marie Mirebeau providing not
only times, but also rays and sensitivity kernels.

The code cand be found at the following address

https://github.com/Mirebeau/HamiltonFastMarching (alias HFM)
============================================================

and illustrations are provided in the context of adaptative grid discretization

https://github.com/Mirebeau/AdaptiveGridDiscretizations

with different notebooks which go far beyond what is needed in seismology.

I recommend the Notebooks_FMM if you are familiar with jupyter tool of training.

The screening of the ADG strategy is beyond this training (but curiosity may
drive you in this automatic differential strategy). The library HFM has been
developped by explicit formulae to be solved numerically.

Therefore, you can download the library HFM (or use conda installed library) and
you can compile it using cmake, following instructions inside this package. 

If you succeed to compile the library, congratulations. Code binaries will be
sitting inside the directory "bin", especially the FileHFM_Isotropic2 (or
FileHFM_Isotropic2.exe under windows)

Such executable code (always rename it FileHFM_Isotropic2.exe if the extension is
missing even under linux) should be copied in the running directory for illustration
how to use such a library through input and output files (agnostic to the computer
langage you are using). I am using it embedded in a reverse-communication (RC)
interface with my Fortan programs for tomography.

For illustration, Python3 codes are provided for showing how to use this RC
interface with the multiplex input and output files.

This strategy is minimal, although I/O issues can be avoided.

Remark:
======

This library has been interfaced with Python and other langages. The Python binding
has to be installed. Direct interface to C++ from Python3 can be from anaconda by
install python and pybind11 (and also numpy and matplotlib)

Such a binding is left to future investigation









