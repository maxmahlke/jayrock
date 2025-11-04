import jayrock
import rocks

bft = rocks.load_bft()
tnos_with_D = bft[
    (bft["sso_class"].str.startswith("KBO")) & (~bft["diameter.value"].isnull())
]
print(
    tnos_with_D[
        ["sso_number", "sso_name", "sso_class", "diameter.value", "albedo.value"]
    ].sort_values("diameter.value", ascending=False)
)
breakpoint()

TARGETS = ["logos", "borasisi-pabu", "sila-nunam", "arrokoth"]

for target in TARGETS:
    # Create target
    tno = jayrock.Target(target)
    tno.compute_ephemeris(cycle=6)

    # Define instrument
    inst = jayrock.Instrument("NIRSpec", mode="IFU")
    inst.disperser = "g395m"
    inst.filter = "f290lp"

    # Define observation date (near vmag min)
    date_obs = tno.get_date_obs(at="vmag_min")

    # Set SNR target
    inst.set_snr_target(100, tno, date_obs)

    # Run simulation
    obs = jayrock.observe(tno, inst, date_obs)

    # Print summary
    obs.print_summary()

    # Plot SNR
    obs.plot_snr()
