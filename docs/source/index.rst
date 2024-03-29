.. Mag2DPoly documentation master file, created by
   sphinx-quickstart on Tue Dec  8 15:57:35 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



``pyMag2DPoly``'s documentation
################################

Overview
*************

**pyMag2DPoly** *is a Python package conceived for forward magnetic anomaly calculation due to two-dimensional polygonal bodies with uniform arbitrary polarization*. 

The formulations implemented in this package are that of Talwani & Heirtzler (1962, 1964), Won & Bevis (1987) and revised Kravchinsky et al. (2019).

If you use this code for research or else, please cite the related paper:

Alessandro Ghirotto, Andrea Zunino, Egidio Armadillo & Klaus Mosegaard (2021). **Magnetic Anomalies Caused by 2D Polygonal Structures with Uniform Arbitrary Polarization: new insights from analytical/numerical comparison among available algorithm formulations**. *Geophysical Research Letters, 48* (7), e2020GL091732.

The specific procedures for each formulation, the analytical/numerical results derived from their comparison and the rectification made to Kravchinsky et al. (2019) algorithm are describde in detail in the paper above.


Installation
============

To download and install the package use "pip3" or "pip" in the following ways: 

.. code-block:: python

    pip3 install Mag2DPoly

or locally, providing the path to the directory Mag2DPoly

.. code-block:: python

    pip3 install path/to/Mag2DPoly

or

.. code-block:: python

    pip3 install https://github.com/inverseproblem/pyMag2DPoly



.. toctree::
   :caption: Contents:
   :maxdepth: 4
   
   userguide.rst
   apilist.rst
   

Indices and tables
********************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
