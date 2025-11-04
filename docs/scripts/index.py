import pickle
from pathlib import Path

# ------
import jayrock

# ------
achilles = jayrock.Target("Achilles")
achilles.compute_ephemeris(cycle=6)

# ------
achilles.print_ephemeris()

# ------
achilles.plot_spectrum(date_obs="2028-03-06")  # specific date

date_thermal_max = achilles.get_date_obs(at="thermal_max")
date_thermal_min = achilles.get_date_obs(at="thermal_min")
achilles.plot_spectrum(
    date_obs=[date_thermal_max, date_thermal_min]
)  # compare at thermal extrema

# ------
# Simulate observation of Achilles on date of minimum Vmag
date_obs = achilles.get_date_obs(at="vmag_min")

observations = []  # list to store observations

inst = jayrock.Instrument("NIRSpec", mode="IFU")

PATH_CACHE = Path("/tmp/obs_nirspec.pkl")

if not PATH_CACHE.is_file():
    for disp, filt in [("G235M", "F170LP"), ("G395M", "F290LP")]:
        # Set filter and disperser
        inst.disperser = disp
        inst.filter = filt

        # Set ngroup, nint, nexp directly
        inst.detector.ngroup = 10  # number of groups per integration
        inst.detector.nint = 1  # number of integrations per exposure
        inst.detector.nexp = 4  # number of exposures -> number of dithers

        obs = jayrock.observe(achilles, inst, date_obs=date_obs)
        observations.append(obs)

    with open(PATH_CACHE, "wb") as f:
        pickle.dump(observations, f)
else:
    with open(PATH_CACHE, "rb") as f:
        observations = pickle.load(f)

# Plot SNR over wavelength
jayrock.plot_snr(observations)

# ------
# # Simulate observation of Achilles on date of minimum Vmag
date_obs = achilles.get_date_obs(at="vmag_min")

observations = []  # list to store observations

inst = jayrock.Instrument("MIRI", mode="MRS")

PATH_CACHE = Path("/tmp/obs_miri_mrs.pkl")

if not PATH_CACHE.is_file():
    # Each channel and disperser is simulated separately
    for channel in ["ch1", "ch2", "ch3", "ch4"]:
        for disp in ["short", "medium", "long"]:
            # Set channel and disperser
            inst.aperture = channel
            inst.disperser = disp

            inst.detector.nexp = 4  # 4-pt dither

            # ------
            # Define SNR target based on channel
            # mrsshort: ch1|ch2 - mrslong: ch3|ch4
            if channel == "ch1":  # 4.9-7.65μm
                snr_target = 10

            if channel == "ch2":  # 7.51-11.70μm
                snr_target = 300

            if channel == "ch3":  # 11.55-17.98μm
                snr_target = 200

            if channel == "ch4":  # 17.70-27.9μm
                snr_target = 100

            # jayrock sets ngroup and nint based on snr target
            inst.set_snr_target(snr_target, achilles, date_obs)

            obs = jayrock.observe(achilles, inst, date_obs=date_obs)
            observations.append(obs)
        with open(PATH_CACHE, "wb") as f:
            pickle.dump(observations, f)
else:
    with open(PATH_CACHE, "rb") as f:
        observations = pickle.load(f)

# Plot SNR over wavelength
jayrock.plot_snr(observations)
