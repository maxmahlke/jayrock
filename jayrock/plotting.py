import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "nirspec": {
        "g140m": {"f100lp": "#0b132b", "f070lp": "#1c2541", "f170lp": "#3a506b"},
        "g235m": {"f170lp": "#092327", "f140m": "#0b5351", "f290lp": "#00a9a5"},
        "g395m": {"f290lp": "#4f6367", "f170lp": "#7a9e9f", "f100lp": "#b8d8d8"},
    },
    "miri": {
        "ch1": {"short": "#0b132b", "medium": "#1c2541", "long": "#3a506b"},
        "ch2": {"short": "#092327", "medium": "#0b5351", "long": "#00a9a5"},
        "ch3": {"short": "#4f6367", "medium": "#7a9e9f", "long": "#b8d8d8"},
        "ch4": {"short": "#064789", "medium": "#427aa1", "long": "#8d99ae"},
    },
}


def get_colors(N, cmap="turbo"):
    """
    Get a list of unique colors.

    Parameters
    ----------
    N : int
        The number of unique colors to return.
    cmap : str
        The matplotlib colormap to sample. Default is 'turbo'

    Returns
    -------
    list of str
        A list of color-hexcodes.

    """
    COLORS = plt.get_cmap(cmap, N)
    return [mpl.colors.rgb2hex(COLORS(i)[:3]) for i in range(N)]


def plot_snr(obs):
    """Plot SNR of observation.

    Parameters
    ----------
    obs : Observation or list of Observation
        Observation object containing target, instrument, and visibility information.
    """

    if not isinstance(obs, list):
        observations = [obs]
    else:
        observations = obs

    fig, ax = plt.subplots()

    for obs in observations:
        inst = obs.instrument.instrument

        # TODO: Simplify this access
        if inst == "miri":
            elem1 = obs.instrument.config["instrument"]["aperture"]
            elem2 = obs.instrument.config["instrument"]["disperser"]
        elif inst == "nirspec":
            elem1 = obs.instrument.config["instrument"]["disperser"]
            elem2 = obs.instrument.config["instrument"]["filter"]

        color = COLORS[inst][elem1][elem2]

        # Indicate element name
        ax.axvline(min(obs.wave), ls="--", color=color, lw=0.5)
        ax.text(
            min(obs.wave) - 0.03,
            0.95 if not obs.faint else 0.05,
            f"{elem1} {elem2} - ngroups={obs.instrument.config['detector']['ngroup']}",
            color=color if not obs.faint else "gray",
            va="top" if not obs.faint else "bottom",
            transform=ax.get_xaxis_transform(),
            rotation=270,
            ha="right",
        )
        # Plot SNR
        ax.scatter(obs.wave, obs.snr, s=3, color=color if not obs.faint else "gray")

        # highlight partially and fully saturated
        ax.scatter(
            obs.wave[obs.partial], obs.snr[obs.partial], s=3, ec="orange", fc="none"
        )
        ax.scatter(obs.wave[obs.full], obs.snr[obs.full], s=3, ec="red", fc="none")

        if np.any(obs.partial):
            ax.text(
                min(obs.wave),
                0.75,
                "partially saturated",
                color="orange",
                va="top",
                transform=ax.get_xaxis_transform(),
                rotation=270,
                ha="right",
            )

        if np.any(obs.full):
            ax.text(
                min(obs.wave),
                0.5,
                "saturated",
                color="red",
                va="top",
                transform=ax.get_xaxis_transform(),
                rotation=270,
                ha="right",
            )

    vmag_min = min(o.vmag for o in observations)
    vmag_max = max(o.vmag for o in observations)
    vmags = (
        "-".join([f"{vmag_min:.2f}", f"{vmag_max:.2f}"])
        if vmag_min != vmag_max
        else f"{vmag_min:.2f}"
    )

    flux_min = min(o.thermal_flux for o in observations)
    flux_max = max(o.thermal_flux for o in observations)
    fluxes = (
        "-".join([f"{flux_min:.2f}", f"{flux_max:.2f}"])
        if flux_min != flux_max
        else f"{flux_min:.2f}"
    )

    ax.set_title(
        f"{inst.upper()} {obs.instrument.mode.upper()} SNR for {obs.target.name} - [{vmags}] Vmag - [{fluxes}] mJy"
    )

    if inst == "miri":
        ax.set_yscale("log")

    ax.set(xlabel=r"Wavelength / \textmu m", ylabel="SNR")
    # fig.savefig(f"{obs.target.name}_{inst}_snr.png", backend="pgf")
    return fig, ax


def compare_observations(obs, which="snr", label=None):
    """Compare observations in a given property."""
    if label is None:
        label = ""

    COLORS = ["black", "orange", "firebrick"]

    fig, ax = plt.subplots()

    for i, ob in enumerate(obs):
        wave, snr = ob.report["1d"]["sn"]

        # Highlight bins with partial/full saturation
        partial = ob.report["1d"]["n_partial_saturated"][1] > 0
        full = ob.report["1d"]["n_full_saturated"][1] > 0

        ax.scatter(wave, snr, color=COLORS[0], s=10)
        ax.scatter(wave[partial], snr[partial], color=COLORS[1], s=10)
        ax.scatter(wave[full], snr[full], color=COLORS[2], s=10)

        ax.text(wave[100], snr[100] + 10, ob.label)

    ax.set(xlabel=r"Wavelength / \textmu m", ylabel="SNR")
    ax.set(title=f"Comparison")
    fig.savefig(ob.PATH_INST / f"compare_{which}_{label}.png")


def plot_spectrum(target, reflected=True, thermal=True, date_obs=None, at=None):
    """Plot the spectrum of the target.

    Parameters
    ----------
    target : Target
        Target to plot the spectrum for.
    reflected : bool, optional
        Whether to plot the reflectance spectrum. Default is True.
    thermal : bool, optional
        Whether to plot the thermal spectrum. Default is True.
    date_obs : str, optional
        Date of observation in ISO format (YYYY-MM-DD). If None, uses the
        date of brightest apparent V magnitude. Default is None.
    at : str or list of str, optional
        Conditions at which to plot the spectrum. Choose from ['vmag_min', 'vmag_max',
        'thermal_max', 'thermal_min']. If None, uses date_obs. Default is None.
    """
    plt.style.use("default")

    # Get list of dates to plot
    dates_obs = []

    if date_obs is not None:
        dates_obs.append((date_obs, ""))

    if not isinstance(at, list) and at is not None:
        at = [at]

    if at is not None:
        VALID_AT = ["vmag_min", "vmag_max", "thermal_max", "thermal_min"]

        for at_ in at:
            if at_ not in VALID_AT:
                raise ValueError(
                    f"Invalid value for at: {at_}. Choose from {VALID_AT}."
                )

            prop, cond = at_.split("_")
            prop = "V" if prop == "vmag" else "thermal_flux"

            if cond == "min":
                idx = target.visibility[prop].idxmin()
            elif cond == "max":
                idx = target.visibility[prop].idxmax()

            date_at = target.visibility.date_obs.values[idx]
            dates_obs.append((date_at, at_))

    if date_obs is None and at is None:
        idx = target.visibility.V.idxmin()
        jayrock.logging.info(f"Plotting spectrum at date of brightest V magnitude.")
        dates_obs.append((target.visibility.date_obs.values[idx], "vmag_min"))

    fig, ax = plt.subplots(constrained_layout=True)

    COLORS_REFL = ["steelblue"] + get_colors(len(dates_obs), cmap="winter")
    COLORS_THERM = ["firebrick"] + get_colors(len(dates_obs), cmap="autumn")
    LS = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]

    for i, (date_obs, at_) in enumerate(dates_obs):
        suffix = f"- {date_obs} - [{at_}]" if at_ != "" else f" - {date_obs}"
        ls = LS[i % len(LS)]
        crefl = COLORS_REFL[i % len(COLORS_REFL)]
        ctherm = COLORS_THERM[i % len(COLORS_THERM)]

        if reflected:
            wave_refl, flux_refl = target.compute_reflectance_spectrum(date_obs)
            ax.plot(wave_refl, flux_refl, label=f"Reflected{suffix}", c=crefl, ls=ls)

        if thermal:
            wave_therm, flux_therm = target.compute_thermal_spectrum(date_obs)
            ax.plot(wave_therm, flux_therm, label=f"Thermal{suffix}", c=ctherm, ls=ls)

    ax.set(xlabel="Wavelength / micron", ylabel="Flux / mJy", xlim=(0, 28))
    ax.set_title(f"Spectrum of ({target.rock.number}) {target.name}")
    ax.legend()
    ax.grid()
    plt.show()
