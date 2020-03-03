.. _installation:

:github_url: https://github.com/ShiveshM/GolemFlavor

************
Installation
************

---
Pip
---

GolemFlavor can be installed using ``pip``:

.. code-block:: bash

    pip install git+https://github.com/ShiveshM/GolemFlavor.git

This installs GolemFlavor, along with all the necessary dependencies such as
NumPy and SciPy.

GolemFlavor uses the IceCube software `GolemFit: The HESE
fitter <https://github.com/IceCubeOpenSource/GolemFit>`_ to fit with IceCube
Astrophysical Flavor data. This software is proprietary and so access is
currently limited to IceCube collaborators. A simple Gaussian likelihood can be
used as a substitute for test purposes if this requirement is not found.

------------
Dependencies
------------

GolemFlavor has the following dependencies:

- `Python <https://www.python.org/>`_ >= 2.7 or >= 3.4
- `NumPy <http://www.numpy.org/>`_
- `SciPy <https://www.scipy.org/>`_
- `Six <https://six.readthedocs.io/>`_
- `mpmath <http://mpmath.org/>`_
- `emcee <https://emcee.readthedocs.io/en/stable/>`_
- `PyMultiNest <https://johannesbuchner.github.io/PyMultiNest/>`_
- `tqdm <https://tqdm.github.io/>`_
- `Shapely <https://shapely.readthedocs.io/en/latest/manual.html>`_
- `Matplotlib <https://matplotlib.org/>`_
- `python-ternary <https://github.com/marcharper/python-ternary>`_
- `GetDist <https://getdist.readthedocs.io/en/latest/>`_

You can use ``pip`` to install the above automatically. Note that ``PyMultiNest``
requires the ``MultiNest`` Bayesian inference library, see [the `PyMultiNest
documentation
<https://johannesbuchner.github.io/PyMultiNest/install.html#prerequisites-for-building-the-libraries>`_
for install instructions.

Additional dependencies:

- `GolemFit <https://github.com/IceCubeOpenSource/GolemFit>`_
