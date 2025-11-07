import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


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


def plot_snr(obs, show=True, save_to=None):
    """Plot SNR of observation.

    Parameters
    ----------
    obs : Observation or list of Observation
        Observation object containing target, instrument, and visibility information.
    show : bool, optional
        Show the figure in the interactive matplotlib window. Default is True.
    save_to : str, optional
        Save the figure to the specified path. Default is None, in which case the figure
        is not saved.

    Returns
    -------
    matplotlib.figure.Figure, matplotlib.axes._axes.Axes
        The figure and axis of the created plot.
    """
    # Increase font size
    mpl.rcParams.update({"font.size": 12})
    # Increase tick length
    mpl.rcParams.update({"xtick.major.size": 8, "ytick.major.size": 8})
    mpl.rcParams.update({"xtick.minor.size": 6, "ytick.minor.size": 6})

    if not isinstance(obs, list):
        observations = [obs]
    else:
        observations = obs

    fig, ax = plt.subplots(constrained_layout=True)
    if len(observations) > 2:
        COLORS = get_colors(len(observations), cmap="managua")[::-1]
    else:
        COLORS = ["steelblue", "firebrick"]

    elems_labels = []

    # ------
    # Add SNRs
    # text_settings = {}
    for obs_idx, obs in enumerate(observations):
        inst = obs.instrument
        elems_labels.append(
            (min(obs.wave), (inst.name, inst.mode, inst.primary, inst.secondary))
        )

        # Add information
        ngroup = f"ngroup={obs.instrument.detector.ngroup}"
        nint = f"nint={obs.instrument.detector.nint}"
        nexp = f"nexp={obs.instrument.detector.nexp}"

        texp = obs.texp if obs.texp < 1000 else obs.texp / 60
        texp = f"texp={texp:.1f}{'s' if obs.texp < 1000 else 'min'}"

        desc = " · ".join([ngroup, nint, nexp, texp])

        # Plot SNR
        ax.scatter(obs.wave, obs.snr_1d, s=3, color=COLORS[obs_idx], label=desc)

        if any(obs.partially_saturated):
            ax.scatter(
                obs.wave[obs.partially_saturated],
                obs.snr_1d[obs.partially_saturated],
                color="orange",
                alpha=0.3,
                zorder=-100,
                marker="s",
                s=40,
            )
            ax.text(
                0.05,
                1.01,
                "! Partial Saturation !",
                ha="left",
                va="bottom",
                color="orange",
                transform=ax.transAxes,
            )

        if any(obs.fully_saturated):
            ax.scatter(
                obs.wave[obs.fully_saturated],
                obs.snr_1d[obs.fully_saturated],
                color="purple",
                alpha=0.3,
                zorder=-100,
                marker="s",
                s=40,
            )
            ax.text(
                0.95,
                1.01,
                "! Full Saturation !",
                ha="right",
                color="purple",
                transform=ax.transAxes,
            )

    # ------
    # For all aperture, disperser, filter pairs, indicate their position
    for wave, (inst, mode, elem1, elem2) in set(elems_labels):
        ax.axvline(wave, ls="-", c="black", lw=0.5)
        opts = dict(va="bottom", rotation=270, ha="left", color="black")
        label = f"{elem1}/{elem2}"
        ax.text(wave, 0.03, label, transform=ax.get_xaxis_transform(), **opts)

    if any(obs.instrument.instrument == "miri" for obs in observations):
        ax.set_yscale("log")
    else:
        ax.legend()  # for non-MIRI observations, show legend. for MIRI, too crowded

    ax.set(xlabel=r"Wavelength / micron", ylabel="SNR")
    ax.grid(which="major", ls="--", c="gray", axis="y", lw=0.4)
    ax.grid(which="minor", ls=":", c="gray", axis="y", lw=0.4)

    if show:
        plt.show()
    if save_to is not None:
        fig.savefig(save_to)
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


def plot_spectrum(
    target, date_obs, reflected=True, thermal=True, show=True, save_to=None
):
    """Plot the spectrum of the target.

    Parameters
    ----------
    target : Target
        Target to plot the spectrum for.
    date_obs : str or list of str
        Date of observation in ISO format (YYYY-MM-DD).
    reflected : bool, optional
        Whether to plot the reflectance spectrum. Default is True.
    thermal : bool, optional
        Whether to plot the thermal spectrum. Default is True.
    show : bool, optional
        Whether to show the plot. Default is True.
    save_to : str, optional
        Path to save the plot. Default is None.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object.
    matplotlib.axes.Axes
        The axes object.
    """
    if not isinstance(date_obs, list):
        date_obs = [date_obs]

    fig, ax = plt.subplots(constrained_layout=True)

    COLORS_REFL = ["steelblue"] + get_colors(len(date_obs), cmap="winter")
    COLORS_THERM = ["firebrick"] + get_colors(len(date_obs), cmap="autumn")
    LS = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]

    for i, date in enumerate(date_obs):
        suffix = f" - {date}"
        ls = LS[i % len(LS)]
        crefl = COLORS_REFL[i % len(COLORS_REFL)]
        ctherm = COLORS_THERM[i % len(COLORS_THERM)]

        if reflected:
            wave_refl, flux_refl = target.compute_reflectance_spectrum(date)
            ax.plot(wave_refl, flux_refl, label=f"Reflected{suffix}", c=crefl, ls=ls)

        if thermal:
            wave_therm, flux_therm = target.compute_thermal_spectrum(date)
            ax.plot(wave_therm, flux_therm, label=f"Thermal {suffix}", c=ctherm, ls=ls)

        ax.plot(wave_therm, flux_therm + flux_refl, c="black", ls=ls, lw=0.7)

    ax.set(xlabel="Wavelength / micron", ylabel="Flux / mJy", xlim=(0, 28))
    ax.set_title(f"Spectrum of ({target.rock.number}) {target.name}")
    if reflected and thermal:
        ax.set_yscale("log")
        ax.set_ylim(bottom=min(flux_refl) * 0.8)
    ax.legend()
    if save_to is not None:
        fig.savefig(save_to)
    if show:
        plt.show()
    return fig, ax
