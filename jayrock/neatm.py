"""Near-Earth Asteroid Thermal Model from Harris (1998).

This code is essentially a copy of the implementation written by Michael S. P. Kelley.
    https://github.com/mkelley/mskpy
"""

from astropy import units as u
import numpy as np
from scipy.integrate import quad


class NEATM:
    """The Near Earth Asteroid Thermal Model."""

    def __init__(self, diameter, albedo, beaming=1.0, emissivity=0.9):
        """Initialize the NEATM model.

        Parameters
        ----------
        diameter : float
          Asteroid diameter in km.
        albedo : float
          Asteroid geometric albedo.
        beaming : float
          Beaming parameter. Optional, default is 1.0.
        emissivity : float
          Emissivity. Optional, default is 0.9.
        """

        self.diameter = diameter
        self.albedo = albedo
        self.beaming = beaming
        self.emissivity = emissivity

    @property
    def A(self):
        """The bond albedo (albedo * phase integral)."""
        G = 0.15  # H,G slope parameter
        A = self.albedo * (0.290 + 0.684 * G)
        return A

    def Tss(self, rh):
        """The sub-solar point temperature.

        Parameters
        ----------
        rh : float
          Heliocentric distance in au.

        Returns
        -------
        float
          The sub-solar temperature in K.
        """

        Fsun = 1367.567 / rh**2  # W / m2
        sigma = 5.670373e-08  # W / (K4 m2)
        Tss = (
            ((1.0 - self.A) * Fsun) / abs(self.beaming) / self.emissivity / sigma
        ) ** 0.25
        return Tss

    def _point_emission(self, phi, theta, wave, Tss):
        """The emission from a single point.

        Parameters
        ----------
        phi : float
            The longitude in radians.
        theta : float
            The latitude in radians.
        wave : float
            The wavelength in microns.
        Tss : float
            The sub-solar temperature in K.

        Returns
        -------
        B : float
            The emission from the point.
        """
        T = Tss * np.cos(phi) ** 0.25 * np.cos(theta) ** 0.25
        B = planck(wave, T)  # W / (m2 sr Hz)
        return B * np.pi * np.cos(phi) ** 2  # W / (m2 Hz)

    def _latitude_emission(self, theta, wave, Tss, phase):
        """The emission from a single latitude.

        Parameters
        ----------
        theta : float
            The latitude in radians.
        wave : float
            The wavelength in microns.
        Tss : float
            The sub-solar temperature in K.
        phase : float
            The phase angle in radians.

        Returns
        -------
        fluxd : float
            The flux density from the latitude band.
        """
        integral = quad(
            self._point_emission,
            0.0,
            np.pi / 2.0,
            args=(theta, wave, Tss),
            epsrel=1e-4,
        )
        fluxd = integral[0] * np.cos(theta - phase)
        return fluxd if not np.isnan(fluxd) else 0.0

    def fluxd(self, wave, geometry):
        """Flux density.

        Parameters
        ----------
        wave : float or list of float
          The wavelengths at which to compute the emission in microns.
        geometry : dict
          A dictionary containing
            rh : heliocentric distance in au
            delta : observer-target distance in au
            phase : phase angle in degrees

        Returns
        -------
        fluxd : Quantity
          The flux density from the whole asteroid in mJy.
        """
        phase = np.radians(geometry["phase"])

        if not isinstance(wave, (list, np.ndarray)):
            wave = np.array([wave])

        Tss = self.Tss(geometry["rh"])
        fluxd = np.zeros(len(wave))

        for i, w in enumerate(wave):
            fluxd[i] = quad(
                self._latitude_emission,
                -np.pi / 2.0 + np.abs(phase),
                np.pi / 2.0,
                args=(w, Tss, np.abs(phase)),
                epsrel=1e-3,
            )[0]

        delta_km = geometry["delta"] * u.au.to(u.km)

        fluxd *= (
            self.emissivity * (self.diameter / (delta_km)) ** 2 / np.pi / 2.0
        )  # W/m^2/Hz

        fluxd = fluxd * u.Unit("W / (m2 Hz)")
        equiv = u.spectral_density(wave * u.um)
        fluxd = fluxd.to("mJy", equivalencies=equiv).value
        if len(fluxd) == 1:
            return fluxd[0]
        else:
            return fluxd


def planck(wave, T):
    """The Planck function.

    Parameters
    ----------
    wave : float or list of float
      The wavelength(s) to evaluate the Planck function in microns.
    T : float or list of float
      The temperature(s) of the Planck function in K.

    Returns
    -------
    B : float
      The Planck function evaluated at `wave` and `T` in units of W/(m2 Hz sr).
    """
    wave = wave * 1e-6  # convert microns to meters

    c1 = 3.9728913665386057e-25  # J m
    c2 = 0.0143877695998  # K m

    # Clip the argument to a safe maximum value (e.g., 700.0)
    # 700.0 is safe from overflow (limit is ~709.78)
    exp_arg = c2 / (wave * T)
    exp_arg = np.clip(exp_arg, a_min=None, a_max=700.0)

    a = np.exp(exp_arg)
    B = c1 / (wave**3 * (a - 1.0))
    return B
