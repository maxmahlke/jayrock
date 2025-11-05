<p align="center">
  <img width="260" src="https://raw.githubusercontent.com/maxmahlke/jayrock/master/docs/_static/logo_jayrock.svg">
</p>

<p align="center">
  <a href="https://github.com/maxmahlke/jayrock#features"> Features </a> - <a href="https://github.com/maxmahlke/jayrock#install"> Install </a> - <a href="https://github.com/maxmahlke/jayrock#documentation"> Documentation </a>
</p>

<div align="center">
  <a href="https://img.shields.io/pypi/pyversions/jayrock">
    <img src="https://img.shields.io/pypi/pyversions/jayrock"/>
  </a>
  <a href="https://img.shields.io/pypi/v/jayrock">
    <img src="https://img.shields.io/pypi/v/jayrock"/>
  </a>
  <a href="https://readthedocs.org/projects/jayrock/badge/?version=latest">
    <img src="https://readthedocs.org/projects/jayrock/badge/?version=latest"/>
  </a>
</div>

JWST observation planning and simulation for minor bodies across the Solar
System. You define the target and instrument(s) – `jayrock` does the rest.

```python
>>> import jayrock
```

## Features

### Computation of target ephemeris

Define a solar system target by name and then compute and display its
JWST-observable ephemeris for a specified cycle, including visibility windows,
magnitude changes, and thermal flux estimates.

``` python
>>> achilles = jayrock.Target("achilles") # define target by name or designation
>>> achilles.compute_ephemeris(cycle=6)   # compute ephemeris for JWST Cycle 6
>>> achilles.print_ephemeris()
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
```

### Modelling of reflectance spectrum and thermal emission

Identify the dates of maximum and minimum thermal emission within the
visibility window and then plot the target's modelled reflectance and thermal
spectrum for a direct comparison between the two observation dates.

```python
>>> # compare spectra at dates of max/min thermal emission
>>> date_thermal_max = achilles.get_date_obs(at="thermal_max")
>>> date_thermal_min = achilles.get_date_obs(at="thermal_min")
>>> achilles.plot_spectrum(date_obs=[date_thermal_max, date_thermal_min])
```

![Modelled spectrum of (588) Achilles at two different observation dates](https://jayrock.readthedocs.io/en/latest/_images/achilles_spec_thermal_dark.png)
### Compute instrument configuration to reach SNR target

Configure a JWST instrument and calculate the required exposure time parameters
(groups and integrations) to achieve a specified Signal-to-Noise Ratio (SNR) at
a particular wavelength.

```python
>>> miri_lrs = jayrock.Instrument('MIRI', mode='LRSSLIT')
>>> miri_lrs.detector.nexp = 2  # 2-pt dither
>>> miri_lrs.detector.readout_pattern = 'fastr1'  # recommended for LRS
>>> # Compute exposure times for SNR target of 300 at 10 micron
>>> miri_lrs.set_snr_target(300, target, date_thermal_min, wave=10)
>>> # -> jayrock defines ngroup and nint for you
```

### Visualise SNR across detector configurations and possible saturation issues

Simulate multiple observations using manually set exposure parameters across
different detector configurations (apertures/filters/dispersers). Plot
the resulting simulated SNR over the full wavelength range for all defined
observations.

```python
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
```


![Simulated SNR of NIRSpec observation of (588) Achilles](https://jayrock.readthedocs.io/en/latest/_images/achilles_nirspec_dark.png)

## Install

Install from PyPi using `pip`:

     $ pip install jayrock

The minimum required `python` version is 3.11. `jayrock` relies on `pandeia`, the `python` engine behind the JWST ETC.
Detailed install instructions for this dependency can be found [here](https://jayrock.readthedocs.io/en/latest/getting_started.html).


## Documentation

Check out the documentation at [jayrock.readthedocs.io](https://jayrock.readthedocs.io/en/latest/).
