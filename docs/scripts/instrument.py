import jayrock

inst = jayrock.Instrument("NIRSpec", mode="IFU")
# inst.print_config()

# print(inst.disperser)
inst.disperser = "g235m"
inst.filter = "f170lp"
# print(inst.disperser)

# inst.detector.ngroup = 5  # number of groups per integration
# inst.detector.nint = 1  # number of integrations per exposure
inst.detector.nexp = 4  # number of exposures -> number of dithers
inst.detector.readout_pattern = "nrsirs2rapid"

# inst.detector.ngroup = 120  # number of groups per integration
# inst.detector.readout_pattern = "nrs"  # number of groups per integration

# Get target and observation date
target = jayrock.Target("Pawlowia")
target.compute_ephemeris(cycle=6)
date_obs = target.get_date_obs(at="vmag_min")

# Using the Instrument config from above (NIRSpec IFU G395M/F290LP)
# SNR = 100  # target SNR to reach
# inst.set_snr_target(SNR, target, date_obs)
# #
# inst.set_snr_target(SNR, target, date_obs, wave=2.7)

SNR_TARGET_RANGE = (80, 120)
inst.set_snr_target(SNR_TARGET_RANGE, target, date_obs, wave=2.7)
