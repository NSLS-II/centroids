Theory of Operation
===================

Contents
--------

.. toctree::

   first_pass
   second_pass

Introduction
------------

The centroids, photon finding routines work through the following principles.
Initially, the code makes a first pass through each image pixel by pixel
looking for potential photon events. This search algorithm is described in
:doc:`first_pass`. Once each image has been searched, and a list of potential
events found, the code then processes each potential photon in turn.

The second pass, described in :doc:`second_pass` evaluates each potential
photon event applying a number of filters to determine if the event is a true
photon. If the event is determined to be a photon, then the code performs
calculations on the photon, storing the values in a final list.

:doc:`first_pass`

:doc:`second_pass`