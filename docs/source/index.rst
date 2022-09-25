.. image:: ../image/banner.png
   :width: 500px
   :height: 500px
   :scale: 50 %
   :alt: alternate text
   :align: left

.. image:: https://readthedocs.org/projects/damascus/badge/?version=latest
	:target: http://damascus.readthedocs.io/en/latest/?badge=latest
	:alt: Documentation Status
.. image:: https://travis-ci.com/temken/DaMaSCUS.svg?branch=master
    :target: https://travis-ci.com/temken/DaMaSCUS
.. image:: https://codecov.io/gh/temken/DaMaSCUS/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/temken/DaMaSCUS
	
==================================================================
DaMaSCUS - Dark Matter Simulation Code for Underground Simulations
==================================================================

Authors: `Timon Emken <https://timonemken.com>`_ & `Chris Kouvaris <http://cp3-origins.dk/people/staff/kouvaris>`_

.. raw:: html

	<a href="http://ascl.net/1706.003"><img src="https://img.shields.io/badge/ascl-1706.003-blue.svg?colorB=262255" alt="ascl:1706.003" /></a>
   
|

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3726878.svg
   :target: https://doi.org/10.5281/zenodo.3726878

* DaMaSCUS is a MC simulator of dark matter particles as they move through the Earth and scatter on terrestrial nuclei. 
* It allows to compute the local distortions of the DM density and velocity distribution caused by collisions with nuclei. 
* The distorted distribution functions and redistributed densities are used to give precise estimates of time-dependent signal rates for direct detection experiments and diurnal modulations.
* A full, realistic model of the Earth is implemented as well as the Earth's time-dependent velocity and orientation in the galactic frame.
* DaMaSCUS is written in C++ and fully parallelized (MPI).

For the underlying physics check out the `paper <https://arxiv.org/abs/1706.02249>`_ . For the code visit the `repository <https://github.com/temken/DaMaSCUS>`_ .


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   gettingstarted
   usage
   citation
   releases
   license
   contact

=============
Visualization
=============
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Simulation of a single DM particle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

	<video controls src="_static/Damascus_Single_Track.mp4"></video> 

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Simulation of a multiple particles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

	<video controls src="_static/Damascus_Multi_Particles.mp4"></video> 

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Evolution of Isodetection Rings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: ../image/IsoRings.gif
   :alt: alternate text
   :align: left

=================
Contact & Support
=================

For questions, bug reports or other suggestions please contact

	`Timon Emken <https://timonemken.com>`_ (emken@chalmers.se)


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`

.. * :ref:`search`
