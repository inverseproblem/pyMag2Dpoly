API, list of functions
########################################

Data structures
===============

.. automodule:: mag2dpoly.magdatastruct 
   :members: BodySegments2D, MagPolyBodies2D, MagnetizVector 
   :imported-members: 

.. warning::
   
   Vertices of the polygonal bodies must be provided 
   counterclockwise to the function ``BodySegments2D``
   to perform magnetic anomaly calculation using the
   functions in the next section **Forward functions**


Forward functions
=================

Single polygonal body
^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mag2dpoly.mag2dpolybodies 
   :members: tmagpoly2D, tmagpoly2Dgen
   :imported-members:
   :noindex:

Multiple polygonal bodies
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mag2dpoly.mag2dpolybodies 
   :members: tmagpolybodies2D, tmagpolybodies2Dgen
   :imported-members:
   :noindex:

Forward algorithms alone
^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mag2dpoly.mag2dpolybodies 
   :members: tmagtalwani, tmagtalwanired, tmagkrav, tmagwonbev
   :imported-members:
   :noindex:

.. note::
   
   These functions are not exported. To call them 
   type ``mag2dpoly.`` before the name of the functions.

 
Useful functions
^^^^^^^^^^^^^^^^

.. automodule:: mag2dpoly.magutils
   :members:
   :imported-members:
   
.. automodule:: mag2dpoly.mag2dpolybodies 
   :members: checkanticlockwiseorder
   :imported-members:
   :noindex:
   
.. note::

   These functions are not exported. To call them
   type ``mag2dpoly.`` before the name of the functions.

