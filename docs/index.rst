###########
``jayrock``
###########

.. raw:: html

    <style> .gray {color:#979eab} </style>

.. role:: gray

JWST observation planning and simulation for minor bodies across the Solar System. You define the target and instrument(s) -- ``jayrock`` does the
rest.\ [#f1]_

.. warning::

   ``jayrock`` relies on existing services and tools, primarily provided directly
   by the STScI. Nevertheless, bugs happen. I strongly recommend to verify the
   planning results obtained here with the tools provided by the STScI (JWST GTVT, ETC, APT).

.. |br| raw:: html

     <br>

.. highlight:: python

This page gives a brief overview.
Refer to the **Planning an Observation** guide on the left for more information and to the **Examples**
for complete workflows.

.. code-block:: python

   import jayrock

:octicon:`telescope;1em` **Compute visibility windows and ephemeris.**

``jayrock`` verifies the target's visibility with `jwst_gtvt <https://github.com/spacetelescope/jwst_gtvt>`_ and retrieves ephemeris
data from `JPL Horizons <https://ssd.jpl.nasa.gov/horizons/>`_.

.. code-block:: python

   achilles = jayrock.Target("achilles") # define target by name or designation
   achilles.compute_ephemeris(cycle=6)   # compute ephemeris for JWST Cycle 6

Get a quick overview of visibility windows with the ``print_ephemeris()`` method.

.. code-block:: python

   achilles.print_ephemeris()

.. dropdown:: Show Output

    .. code-block:: python

      (588) Achilles: Ephemeris from 2027-07-01 to 2028-06-30
     ├── Window 1: 2027-07-02 -> 2027-08-10
     │   ├── Duration         40 days
     │   ├── Vmag             16.18 -> 16.48
     │   └── Thermal @ 15um   393.39 -> 308.31 mJy
     ├── Window 2: 2028-03-02 -> 2028-04-21
     │   ├── Duration         51 days
     │   ├── Vmag             16.54 -> 16.15
     │   └── Thermal @ 15um   279.80 -> 360.74 mJy
     └── errRA/errDec in arcsec: 0.027 / 0.005

|br|
:octicon:`beaker;1em` **Model the target's spectrum at different dates of observation.**

The spectrum is modelled as the sum of reflected sunlight and thermal emission using ``synphot`` and
NEATM. Best-estimate parameters are retrieved under the hood with ``rocks``. To visualise the
spectrum at different dates, use the ``plot_spectrum()`` method.

.. code-block:: python

   # plot spectrum on a specific date
   achilles.plot_spectrum(date_obs="2028-03-06")

   # compare spectra at dates of max/min thermal emission
   date_thermal_max = achilles.get_date_obs(at="thermal_max")
   date_thermal_min = achilles.get_date_obs(at="thermal_min")
   achilles.plot_spectrum(date_obs=[date_thermal_max, date_thermal_min])

.. dropdown:: Show Output

  *On March 6, 2028:*

  .. image:: gfx/achilles_spec_date.png
     :class: only-light
     :align: center
     :width: 60%

  .. image:: gfx/achilles_spec_date_dark.png
     :class: only-dark
     :align: center
     :width: 60%

  |br|

  *At the lowest and highest predicted thermal emission:*

  .. image:: gfx/achilles_spec_thermal.png
     :class: only-light
     :align: center
     :width: 60%

  .. image:: gfx/achilles_spec_thermal_dark.png
     :class: only-dark
     :align: center
     :width: 60%

|br|
:octicon:`zap;1em` **Compute SNRs for many instrument configurations. Define SNR targets.**

Simulate observations with different instruments, modes, filters and
dispersers, conveniently in a loop. Specify exposure times or define SNR
target values and let ``jayrock`` configure the detector. Computations run
locally using ``pandeia`` -- no need for the browser.

.. tab-set::

   .. tab-item:: Defining exposure times

     In this example, we define exposure times directly via the ``ngroup``, ``nint``, and ``nexp``
     settings.

     .. code-block:: python

        # Simulate observation of Achilles on date of minimum Vmag
        date_obs = achilles.get_date_obs(at="vmag_min")

        observations = []  # list to store observations

        inst = jayrock.Instrument("NIRSpec", mode="IFU")

        for disp, filt in [("G235M", "F170LP"), ("G395M", "F290LP")]:

            # Set filter and disperser
            inst.disperser = disp
            inst.filter = filt

            # Set ngroup, nint, nexp directly
            inst.detector.ngroup = 10  # number of groups per integration
            inst.detector.nint = 1     # number of integrations per exposure
            inst.detector.nexp = 4     # number of exposures -> number of dithers

            obs = jayrock.observe(achilles, inst, date_obs=date_obs)
            observations.append(obs)

        # Plot SNR over wavelength
        jayrock.plot_snr(observations)

     .. dropdown:: Show Figure

       *Simulated SNR over wavelength for NIRSpec IFU. Click to zoom.*

       .. image:: gfx/achilles_nirspec.png
          :class: only-light
          :align: center
          :width: 100%

       .. image:: gfx/achilles_nirspec_dark.png
          :class: only-dark
          :align: center
          :width: 100%

   .. tab-item:: Setting SNR targets

      In this example, we define SNR targets instead of exposure times.\ [#f2]_

      .. code-block:: python

         # Simulate observation of Achilles on date of minimum Vmag
         date_obs = achilles.get_date_obs(at="vmag_min")

         observations = []  # list to store observations

         inst = jayrock.Instrument("MIRI", mode="MRS")

         # Each aperture and disperser is simulated separately
         for aperture, disp in zip(['ch1', 'ch2', 'ch3', 'ch4'], ['short', 'medium', 'long']):

             # Set aperture and disperser
             inst.aperture = aperture
             inst.disperser = disp

             inst.detector.nexp = 4  # 4-pt dither

             # ------
             # Define SNR target based on aperture
             # mrsshort: ch1|ch2 - mrslong: ch3|ch4
             if aperture == 'ch1':  # 4.9-7.65μm
               snr_target = 10

             if aperture == 'ch2':  # 7.51-11.70μm
               snr_target = 300

             if aperture == 'ch3':  # 11.55-17.98μm
               snr_target = 200

             if aperture == 'ch4':  # 17.70-27.9μm
               snr_target = 100

             # jayrock sets ngroup and nint based on snr target
             inst.set_snr_target(snr_target, achilles, date_obs)

             obs = jayrock.observe(achilles, inst, date_obs=date_obs)
             observations.append(obs)

         # Plot SNR over wavelength
         jayrock.plot_snr(observations)

      .. dropdown:: Show Figure

        *Simulated SNR over wavelength for MIRI MRS. All channels have the requested SNRs at their reference wavelengths. Click to zoom.*

        .. image:: gfx/achilles_miri.png
           :class: only-light
           :align: center
           :width: 100%

        .. image:: gfx/achilles_miri_dark.png
           :class: only-dark
           :align: center
           :width: 100%


.. rubric:: Footnotes

.. [#f1] Latest version: 0.0.2  - `What's new? <https://github.com/maxmahlke/jayrock/blob/master/CHANGELOG.md>`_  | Bug or feature request? Open an issue on `GitHub <https://github.com/maxmahlke/jayrock/issues>`_ or send me an email.

.. [#f2] This example does not account for the fact that channels 1/2 and 3/4 use the same detectors and therefore need to have the same ``ngroup``/``nint`` settings in their ``short``, ``medium``, and ``long`` dispersers.

.. toctree::
   :maxdepth: 2
   :caption: Welcome to jayrock
   :hidden:

   Home<self>
   Installation<getting_started>

.. toctree::
   :maxdepth: 2
   :caption: Planning an Observation
   :hidden:

   Identify your Target<target>
   Configure the Instrument<instrument>
   Simulate the Observation<observation>

.. toctree::
   :maxdepth: 2
   :caption: Examples
   :hidden:

   Vesta family with MIRI<vesta_miri>
   TNOs with NIRSpec<tno_nirspec>

.. toctree::
   :maxdepth: 2
   :caption: Reference
   :hidden:

   API Reference<api>
