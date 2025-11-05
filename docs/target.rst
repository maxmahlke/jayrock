Identify your Target
====================

The ``jayrock.Target`` class handles identification, gathering physical
properties, and calculating orbital visibility/fluxes (ephemeris) for your
minor body as seen from JWST.

Defining a Target
-----------------

Identify your target by passing a name, designation, or number. ``jayrock`` uses `rocks <https://rocks.readthedocs.io/en/latest/cli.html>`_ to automatically retrieve available physical parameters.


.. code-block:: python

   import jayrock
   luisa = jayrock.Target('Luisa')        # (599) Luisa
   chariklo = jayrock.Target('1997 CU26') # (10199) Chariklo
   pluto = jayrock.Target(134340)         # (134340) Pluto


Get ephemeris as seen from JWST
-------------------------------

After defining the target, we want to know whether it is observable from JWST during a given timeframe (cycle)
and at what visual magnitude / thermal flux we can observe it. The ``compute_ephemeris()`` method computes visibility windows and fluxes using `jwst_gtvt <https://github.com/spacetelescope/jwst_gtvt>`_ and `JPL Horizons <https://ssd.jpl.nasa.gov/horizons/>`_.

.. code-block:: python

   # Ephemeris of (599) Luisa for Cycle 6 (2027-07-01 - 2028-06-30)
   luisa.compute_ephemeris(cycle=6)

.. TODO: Include mkpy license in package

The result of the query is stored in the ``ephemeris`` attribute of the target, a `pandas.DataFrame <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html>`_. Each row in the dataframe represents one date on which the target is visible and contains viewing geometry, fluxes, and more.
To get a high-level overview, use the ``print_ephemeris()`` method:

.. code-block:: python

   >>> luisa.print_ephemeris()
   (599) Luisa: Ephemeris from 2027-07-01 to 2028-06-30
   ├── Window 1: 2027-07-12 -> 2027-09-17
   │   ├── Duration         68 days
   │   ├── Vmag             12.72 -> 11.68
   │   └── Thermal @ 15um   12657.03 -> 26584.46 mJy
   ├── Window 2: 2027-12-02 -> 2028-01-31
   │   ├── Duration         61 days
   │   ├── Vmag             12.18 -> 13.55
   │   └── Thermal @ 15um   16169.68 -> 5146.10 mJy
   └── errRA/errDec in arcsec: 0.015 / 0.013


.. warning::

   Visibility windows can slightly vary between those obtained by the APT and
   those obtained here with ``jwst_mtvt`` (order of one day). The APT has the
   final authority.


Select a date of observation
----------------------------

Simulating an observation requires a date of observation ``date_obs``. You can retrieve all possible dates withe the ``get_date_obs()`` method or select one based on brightness criteria using the ``at`` argument.

.. code-block:: python

   >>> dates_obs = luisa.get_date_obs()  # gets all valid date_obs
   >>> print(dates_obs)
   ['2027-07-12', '2027-07-13', '2027-07-14', ..., '2028-01-29', '2028-01-30', '2028-01-31']
   >>> date_obs = luisa.get_date_obs(at="vmag_min")  # gets date_obs at minimum V magnitude
   >>> print(date_obs)
   '20207-09-17'

Valid conditions are ``vmag_min``, ``vmag_max``, ``thermal_min``, ``thermal_max``.

.. note::

   ``vmag_min`` refers to the numerical minimum of V, hence the brightest apparent magnitude.

Modeling the spectrum
---------------------

Target spectra are modelled as superposition of a reflected and a thermal component.

.. _refl:

Reflected Component
^^^^^^^^^^^^^^^^^^^

The reflected component of the minor planet is modeled using a
`Phoenix stellar model
<https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/phoenix-models-available-in-synphot>`_
of a G2V star. The spectrum is normalised to the apparent Johnson V magnitude
as seen from JWST on the date of observation. This setup can
be replicated in the `ETC <https://jwst.etc.stsci.edu/>`_
in the following way:

.. dropdown:: Show Source Editor

    .. tab-set::

       .. tab-item:: Source Editor - Continuum

            .. image:: gfx/etc_source_light.png
               :class: only-light
               :align: center
               :width: 600

            .. image:: gfx/etc_source_dark.png
               :class: only-dark
               :align: center
               :width: 600

       .. tab-item:: Source Editor - Renorm

            .. image:: gfx/etc_norm_light.png
               :class: only-light
               :align: center
               :width: 600

            .. image:: gfx/etc_norm_dark.png
               :class: only-dark
               :align: center
               :width: 600


.. _thermal:

Thermal Component
^^^^^^^^^^^^^^^^^

The thermal component is modelled via the Near-Earth Asteroid Thermal Model (NEATM) by `Harris (1998) <https://ui.adsabs.harvard.edu/abs/1998Icar..131..291H/abstract>`_. The implementation in ``jayrock`` is essentially a copy of the code in Michael Kelley's `mskpy <https://github.com/mkelley/mskpy>`_.
The following parameters are required to model the target's spectrum. If the parameter is unknown (i.e. not available via ``rocks``), ``jayrock`` assumes the given default value.

+-------------------------------+---------------------------------------------------------------------------+-------------+
| **Name**                      | **Description**                                                           | **Default** |
+-------------------------------+---------------------------------------------------------------------------+-------------+
| albedo                        | Target geometrical albedo in V.                                           | 0.1         |
+-------------------------------+---------------------------------------------------------------------------+-------------+
| diameter                      | Target diameter in km.                                                    | 30          |
+-------------------------------+---------------------------------------------------------------------------+-------------+
| beaming (:math:`\eta`)        |   Parameter accounting for surface roughness and thermal inertia effects. | 1.0         |
+-------------------------------+---------------------------------------------------------------------------+-------------+
| emissivity (:math:`\epsilon`) |  Ratio of asteroid's thermal radiation to ideal blackbody's radiation.    | 0.9         |
+-------------------------------+---------------------------------------------------------------------------+-------------+

You can override the default values using the dot-notation.

.. code-block:: python

   luisa.albedo = 0.11
   luisa.diameter = 68
   luisa.beaming = 1
   luisa.emissivity = 1

Exporting and plotting
^^^^^^^^^^^^^^^^^^^^^^

Export the full spectrum or individual components (e.g. for verification with the online ETC):

.. code-block:: python

  luisa.export_spectrum('luisa_spectrum.txt', date_obs='2027-12-24')
  luisa.export_spectrum('luisa_spectrum_refl.txt', date_obs='2027-12-24', thermal=False)  # only reflected component
  luisa.export_spectrum('luisa_spectrum_thermal.txt', date_obs='2027-12-24', reflected=False)  # only thermal component


Plot the spectrum of a target at a given ``date_obs``:

.. code-block:: python

  luisa.plot_spectrum(date_obs='2027-12-24')

.. dropdown:: Show Figure

  .. image:: gfx/luisa_spectrum.png
     :class: only-light
     :align: center
     :width: 600

  .. image:: gfx/luisa_spectrum_dark.png
     :class: only-dark
     :align: center
     :width: 600

You can provide a list of ``date_obs`` to compare the fluxes at different dates.

.. code-block:: python

  date_vmag_min = luisa.get_date_obs(at='vmag_min')
  date_vmag_max = luisa.get_date_obs(at='vmag_max')
  luisa.plot_spectrum(date_obs=[date_vmag_min, date_vmag_max])

.. dropdown:: Show Figure

  .. image:: gfx/luisa_spectra.png
     :class: only-light
     :align: center
     :width: 600

  .. image:: gfx/luisa_spectra_dark.png
     :class: only-dark
     :align: center
     :width: 600
