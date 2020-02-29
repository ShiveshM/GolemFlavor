# GolemFlavor

GolemFlavor is a Python package for running an analysis pipeline using
`GolemFit`.

![GolemFlavor Logo](logo.png)

## Installation
GolemFlavor can be installed using `pip`
```
pip install git+https://github.com/ShiveshM/flavour_ratio.git
```
This installs GolemFlavor, along with all the necessary dependencies such as
NumPy and SciPy.

GolemFlavor uses the IceCube software [`GolemFit: The HESE
fitter`](https://github.com/IceCubeOpenSource/GolemFit) to fit with IceCube
HESE data. Current access is limited to IceCube collaborators. A simple
Gaussian likelihood can be used instead for test purposes if this requirement
is not found.

### Dependencies

GolemFlavor has the following dependencies:
* [`Python`](https://www.python.org/) 2.7
* [`NumPy`](http://www.numpy.org/)
* [`SciPy`](https://www.scipy.org/)
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
documentation]](https://johannesbuchner.github.io/PyMultiNest/install.html#prerequisites-for-building-the-libraries)
for install instructions.

Additional dependencies:
* [`GolemFit`](https://github.com/IceCubeOpenSource/GolemFit)

## License

[MIT License](LICENSE)

Copyright (c) 2020 Shivesh Mandalia
