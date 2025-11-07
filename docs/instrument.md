# Set up your Instrument

The ``jayrock.Instrument`` class allows you to define and configure the JWST
instrument and observing mode used for your observation. This setup determines
the specific optics (filters, gratings) and detector settings (exposure time,
readouts) used for the observation.

## Choose the instrument

To begin, you select the instrument (e.g. NIRSpec) and its observing mode (e.g. IFU).
Refer to the ``pandeia`` [documentation](https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-pandeia-engine-tutorial/pandeia-quickstart#PandeiaQuickstart-Appendices) for supported instruments and modes.


```python
   import jayrock
   inst = jayrock.Instrument("NIRSpec", mode="IFU")
```

Your instrument is configured with default optical and detector settings.
The complete configuration is stored as a nested dictionary in the ``inst.config`` attribute. Use ``inst.print_config()`` to view all settings:


```python
   >>> inst.print_config()
```

```python
   {
       "instrument": {
           "aperture": "ifu",
           "disperser": "g140m",
           "filter": "f100lp",
           "instrument": "nirspec",
           "mode": "ifu"
       },
       "detector": {
           "nexp": 1,
           "ngroup": 10,
           "nint": 1,
           "readout_pattern": "nrs",
           "subarray": "full"
       }
   }
```

## Configure the optics

Once the instrument is defined, you can configure the optical components (aperture, disperser (grating), and filter) using the dot notation.

```python
  inst.aperture = "ifu"
  inst.disperser = "g235m"
  inst.filter = "f170lp"
```

Refer to the ``pandeia`` [documentation](https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-pandeia-engine-tutorial/pandeia-quickstart#PandeiaQuickstart-Appendices) for dispersers and filters per instrument and mode.

## Configure the detector

The detector settings define the number of integrations (= exposure time) and read-out
parameters. You can either define the integration parameters yourself or define
an SNR target and let ``jayrock`` determine the required number of groups and
integrations.

:::{margin}
``jayrock`` compares your detector settings against official JWST
recommendations. If you set a discouraged value, a warning will be issued
to guide you towards better observing strategies.
:::

### Setting exposure parameters

Directly set integration parameters and read-outs via the ``inst.detector`` dictionary.

```python
  inst.detector.ngroup = 5  # number of groups per integration
  inst.detector.nint = 1    # number of integrations per exposure
  inst.detector.nexp = 4    # number of exposures -> number of dithers
  inst.detector.readout_pattern = "nirs2rapid"
```

:::{margin}
Unlike in the ETC, there is no explicit way to define the dithering
strategy. However, treating ``nexp`` as ``ndit`` gives identical results.
:::

### Setting SNR targets

Instead of manually setting ``ngroup`` and ``nint``,
you can define a target Signal-to-Noise Ratio (SNR) for your observation.
``jayrock`` then determines the minimum ``ngroup`` and ``nint`` parameters to
reach this SNR via a binary search.

This estimation uses your ``Target``'s spectrum and the selected instrument configuration.

```python
  >>> SNR_TARGET = 100

  >>> # Get target and observation date
  >>> target = jayrock.Target("Pawlowia")
  >>> target.compute_ephemeris(cycle=6)
  >>> date_obs = target.get_date_obs(at='vmag_min')

  >>> # Using the Instrument configuration from above (NIRSpec IFU G235M/F170LP)
  >>> inst.set_snr_target(SNR_TARGET, target, date_obs)
  INFO     [jayrock] Searching for minimum ngroup|nint to reach SNR=100.00 at 2.42μm
  INFO     [jayrock]   nint=1 | ngroup=100 -> SNR=390.4 | Texp=98.2min
  INFO     [jayrock]   nint=1 | ngroup=52  -> SNR=299.5 | Texp=51.5min
  INFO     [jayrock]   nint=1 | ngroup=28  -> SNR=222.2 | Texp=28.2min
  INFO     [jayrock]   nint=1 | ngroup=16  -> SNR=161.1 | Texp=16.5min
  INFO     [jayrock]   nint=1 | ngroup=10  -> SNR=115.3 | Texp=10.7min
  INFO     [jayrock]   nint=1 | ngroup=7   -> SNR=84.1  | Texp=7.8min
  INFO     [jayrock]   nint=1 | ngroup=9   -> SNR=105.7 | Texp=9.7min
  INFO     [jayrock]   nint=1 | ngroup=8   -> SNR=95.3  | Texp=8.8min
  INFO     [jayrock] Done. Setting ngroup=9, nint=1.
```

``Texp`` indicates the total exposure time of the observation, accounting for ``ngroup``, ``nint``, and ``nexp``.

By default, the SNR is evaluated at the [midpoint wavelength](https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-calculations-page-overview/jwst-etc-wavelength-of-interest#gsc.tab=0) of the selected
aperture/disperser/filter combination (following the online ETC behaviour). You can override this using the ``wave`` argument:

```python
  >>> inst.set_snr_target(SNR_TARGET, target, date_obs, wave=2.7)
  INFO     [jayrock] Searching for minimum ngroup|nint to reach SNR=100.00 at 2.70μm
  INFO     [jayrock]   nint=1 | ngroup=100 -> SNR=328.3 | Texp=98.2min
  INFO     [jayrock]   nint=1 | ngroup=52  -> SNR=248.7 | Texp=51.5min
  INFO     [jayrock]   nint=1 | ngroup=28  -> SNR=181.7 | Texp=28.2min
  INFO     [jayrock]   nint=1 | ngroup=16  -> SNR=128.7 | Texp=16.5min
  INFO     [jayrock]   nint=1 | ngroup=10  -> SNR=89.2  | Texp=10.7min
  INFO     [jayrock]   nint=1 | ngroup=13  -> SNR=110.6 | Texp=13.6min
  INFO     [jayrock]   nint=1 | ngroup=12  -> SNR=103.9 | Texp=12.6min
  INFO     [jayrock]   nint=1 | ngroup=11  -> SNR=96.8  | Texp=11.7min
  INFO     [jayrock] Done. Setting ngroup=12, nint=1.
```

The minimisation process first varies ``nint`` to find a suitable value, and
then runs several simulations varying ``ngroup`` (7 in total assuming
``ngroup`` is between 5 and 100, the default) to find the minimum value that
meets the SNR target. This takes a few minutes. To speed up the process, you
can define an SNR range instead of a fixed value by providing a tuple
``(SNR_min, SNR_max)``. This way, the search stops as soon as the SNR falls
within the defined range.

```python
  >>> SNR_TARGET_RANGE = (80, 120)
  >>> inst.set_snr_target(SNR_TARGET, target, date_obs, wave=2.7)
  INFO     [jayrock] Searching for minimum ngroup|nint to reach SNR range 80.0-120.0 at 2.70μm
  INFO     [jayrock]   nint=1 | ngroup=5   -> SNR=43.2  | Texp=5.8min
  INFO     [jayrock]   nint=1 | ngroup=100 -> SNR=328.3 | Texp=98.2min
  INFO     [jayrock]   nint=1 | ngroup=52  -> SNR=248.7 | Texp=51.5min
  INFO     [jayrock]   nint=1 | ngroup=28  -> SNR=181.7 | Texp=28.2min
  INFO     [jayrock]   nint=1 | ngroup=16  -> SNR=128.7 | Texp=16.5min
  INFO     [jayrock]   nint=1 | ngroup=10  -> SNR=89.2  | Texp=10.7min
  INFO     [jayrock] Done. Setting ngroup=10, nint=1.
```
