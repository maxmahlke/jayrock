import jayrock

luisa = jayrock.Target("Luisa")  # (599) Luisa
chariklo = jayrock.Target("1997 CU26")  # (10199) Chariklo
pluto = jayrock.Target(134340)  # (134340) Pluto

# Ephemeris of (599) Luisa for Cycle 6 (2027-07-01 - 2028-06-30)
luisa.compute_ephemeris(cycle=6)


print(luisa.ephemeris)
luisa.print_ephemeris()

dates_obs = luisa.get_date_obs()
date_obs = luisa.get_date_obs(at="vmag_min")

print(dates_obs)
print(date_obs)

luisa.emissivity = 1
luisa.beaming = 1
luisa.diameter = 68

luisa.export_spectrum("luisa_spectrum.txt", date_obs="2027-09-17")
luisa.export_spectrum(
    "luisa_spectrum_refl.txt", date_obs="2027-12-24", thermal=False
)  # only reflected component
luisa.export_spectrum(
    "luisa_spectrum_thermal.txt", date_obs="2027-12-24", reflected=False
)  # only thermal component
luisa.plot_spectrum(date_obs="2027-12-24")

date_vmag_min = luisa.get_date_obs(at="vmag_min")
date_vmag_max = luisa.get_date_obs(at="vmag_max")
luisa.plot_spectrum(date_obs=[date_vmag_min, date_vmag_max])
