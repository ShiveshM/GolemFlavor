# GolemFlavor

[![Build Status](https://api.travis-ci.org/ShiveshM/GolemFlavor.svg?branch=master)](https://travis-ci.org/ShiveshM/GolemFlavor)
![Python Version](https://img.shields.io/badge/python-2.7+|3.4+-blue.svg)
[![license](https://img.shields.io/github/license/ShiveshM/GolemFlavor 'license')](https://github.com/ShiveshM/GolemFlavor/blob/master/LICENSE)

GolemFlavor is a Python package for running a Bayesian inference analysis
pipeline using Astrophysical Flavor data taken at
[IceCube](https://icecube.wisc.edu/).

![GolemFlavor Logo](logo.png)

## Overview

### What is Astrophysical Flavor data?
This is data of the *flavor* of a
[neutrino](https://icecube.wisc.edu/outreach/neutrinos) taken at the [IceCube
neutrino observatory](https://icecube.wisc.edu/), which is a cubic kilometer
array of optical sensors embedded in the glacial ice at the South Pole. In
particular, *astrophysical* neutrinos are ones that are very-high-energy and
come from astrophysical origins such as [black
holes](https://doi.org/10.1126/science.aat2890). For more on the physics behind
neutrinos see TODO.

### What does the GolemFlavor package do?
This package provides utilities for astrophysical neutrino propagation and
Bayesian statistical modeling focused on advanced Markov Chain Monte Carlo
(MCMC) algorithms. It has been used to make constraints on New Physics models
in the Astrophysical Flavor, as motivated by the paper [*New Physics in
Astrophysical Neutrino
Flavor*](https://doi.org/10.1103/PhysRevLett.115.161303). For more information
on the statistical modeling see TODO.

## Features
* **Portable Flavor Functions**: A set of useful functions for calculating
    measured flavor compositions given a source composition and a mixing
    matrix.
* **MCMC Algorithms**: Affine invariant and nested sampling algorithms provided
    by [emcee](https://emcee.readthedocs.io/) and
    [MultiNest](https://doi.org/10.1111/j.1365-2966.2009.14548.x).
* **Anarchic Sampling**: Sampling of the neutrino mixing matrix is done
    under the [*neutrino mixing
    anarchy*](https://doi.org/10.1016/j.physletb.2003.08.045) hypothesis to
    ensure an unbiased prior.
* **Distributed and parallel computing**: Scripts included to manage the
    workload across a CPU cluster using
    [HTCondor](https://research.cs.wisc.edu/htcondor/).
* **Visualization**: Produce ternary plots of the flavour composition using the
    [python-ternary](https://zenodo.org/badge/latestdoi/19505/marcharper/python-ternary)
    package and joint posterior plots for analyzing MCMC chains using the
    [getdist](https://getdist.readthedocs.io/en/latest/) package.

## Getting Started


## Installation
GolemFlavor can be installed using `pip`
```
pip install git+https://github.com/ShiveshM/GolemFlavor.git
```
This installs GolemFlavor, along with all the necessary dependencies such as
NumPy and SciPy.

GolemFlavor uses the IceCube software [`GolemFit: The HESE
fitter`](https://github.com/IceCubeOpenSource/GolemFit) to fit with IceCube
Astrophysical Flavor data. This software is proprietary and so access is
currently limited to IceCube collaborators. A simple Gaussian likelihood can be
used as a substitute for test purposes if this requirement is not found.

### Dependencies

GolemFlavor has the following dependencies:
* [`Python`](https://www.python.org/) >= 2.7 or >= 3.4
* [`NumPy`](http://www.numpy.org/)
* [`SciPy`](https://www.scipy.org/)
* [`Six`](https://six.readthedocs.io/)
* [`mpmath`](http://mpmath.org/)
* [`emcee`](https://emcee.readthedocs.io/en/stable/)
* [`PyMultiNest`](https://johannesbuchner.github.io/PyMultiNest/)
* [`tqdm`](https://tqdm.github.io/)
* [`Shapely`](https://shapely.readthedocs.io/en/latest/manual.html)
* [`Matplotlib`](https://matplotlib.org/)
* [`python-ternary`](https://github.com/marcharper/python-ternary)
* [`GetDist`](https://getdist.readthedocs.io/en/latest/)

You can use `pip` to install the above automatically. Note that `PyMultiNest`
requires the `MultiNest` Bayesian inference library, see [the `PyMultiNest`
documentation](https://johannesbuchner.github.io/PyMultiNest/install.html#prerequisites-for-building-the-libraries)
for install instructions.

Additional dependencies:
* [`GolemFit`](https://github.com/IceCubeOpenSource/GolemFit)

## License

[MIT License](LICENSE)

Copyright (c) 2020 Shivesh Mandalia
