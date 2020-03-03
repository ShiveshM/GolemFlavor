.. _overview:

:github_url: https://github.com/ShiveshM/GolemFlavor

********
Overview
********

----------------------------------
What is Astrophysical Flavor data?
----------------------------------

This is data of the *flavor* of a neutrino taken at the `IceCube neutrino
observatory <https://icecube.wisc.edu/>`_, which is a cubic kilometer array of
optical sensors embedded in the glacial ice at the South Pole. In particular,
*astrophysical* neutrinos are ones that are very-high-energy and come from
astrophysical origins such as `active galactic nuclei
<https://doi.org/10.1126/science.aat2890>`_. For more on the physics behind
neutrinos see the :doc:`physics` section.

-------------------------------------
What does the GolemFlavor package do?
-------------------------------------

This package provides utilities for astrophysical neutrino propagation and
Bayesian statistical modeling focused on advanced Markov Chain Monte Carlo
(MCMC) algorithms. It has been used to make constraints on New Physics models
in the Astrophysical Flavor, as motivated by the paper `*New Physics in
Astrophysical Neutrino Flavor*
<https://doi.org/10.1103/PhysRevLett.115.161303>`_.  For more information on
the statistical modeling see the :doc:`statistics` section.

--------
Features
--------

- **Portable Flavor Functions**: A set of useful functions for calculating measured flavor compositions given a source composition and a mixing matrix.
- **MCMC Algorithms**: Affine invariant and nested sampling algorithms provided by `emcee <https://emcee.readthedocs.io/>`_ and `MultiNest <https://doi.org/10.1111/j.1365-2966.2009.14548.x>`_.
- **Anarchic Sampling**: Sampling of the neutrino mixing matrix is done under the `*neutrino mixing anarchy* <https://doi.org/10.1016/j.physletb.2003.08.045>`_ hypothesis to ensure an unbiased prior.
- **Distributed and parallel computing**: Scripts included to manage the workload across a CPU cluster using `HTCondor <https://research.cs.wisc.edu/htcondor/>`_.
- **Visualization**: Produce ternary plots of the flavour composition using the `python-ternary <https://zenodo.org/badge/latestdoi/19505/marcharper/python-ternary>`_ package and joint posterior plots for analyzing MCMC chains using the `getdist <https://getdist.readthedocs.io/en/latest/>`_ package.
