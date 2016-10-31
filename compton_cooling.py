import numpy as np

from webbink import radius

MSOL = 2e33 # mass of Sun, g
MU = 0.63 #mean molecular weight
MH = 1.67e-24 # mass of hydrogen, g
G = 6.67e-8
K = 1.38e-16 # Boltzmann Constant, erg/K
KTOKEV = 8.6173e-8 # convert Kelvin to keV
C = 2e10 # speed of light, cm/s
MDOTGS = 6.34e+25


def tshock(mwd, unit='K'):
    '''Returns the shock temperature for the innermost Keplerian orbit of an
    accretion disk around a white dwarf of mass mwd in solar units.
    Units can be Kelvin ('K') or keV ('keV').'''
    # mass in solar masses
    rwd = radius(mwd)
    mwdg = mwd * MSOL
    tshock = 3.0 * MU * MH * G * mwdg / (16.0 * K * rwd)
    if unit.lower() == 'k':
        return tshock
    elif unit.lower() == 'kev':
        return tshock * KTOKEV
    else:
        print("Units must be Kelvin (K) or keV (keV).  Please check your units.")

def urad(lrad, mwd, frad):
    '''Returns the radiation energy density for a radiation source of 
    luminosity Lrad (in erg/s) assuming it covers a fraction frad of
    the surface of a white dwarf of mass mwd.
    Mass should be in solar units.'''
    if frad <= 0 or frad > 1:
        print("frad must be between 0 and 1. Please check and try again.")
    else:
        rwd = radius(mwd)
        urad = lrad / (4.0 * np.pi * rwd ** 2.0 * frad * C)
        return urad

def vkep(mwd):
    rwd = radius(mwd)
    mwdg = mwd * MSOL
    vkep = (G * mwdg / rwd) ** 0.5
    return vkep 

def nebl(mdot, mwd, facc, vacc):
    '''Returns the density in the boundary layer for a given accretion 
    rate mdot, white dwarf mass mwd, accretion velocity vacc (the radial
    velocity of material through the boundary layer) and surface area
    fraction facc over which material accretes on to the white dwarf.
    Mdot should be in solar masses per year, mwd in solar units, and
    vacc in km/s.  Facc should be between >0 and <= 1.'''
    if facc <=0 or facc > 1:
        print("facc should be between >0 and <= 1. Please check and try again.")
    else:
        rwd = radius(mwd)
        mdot = mdot * MDOTGS
        vacc = vacc * 1e5
        nebl = mdot / (4 * np.pi * MU * MH * rwd ** 2 * facc * vacc)
        return nebl

def facc(mwd, mdot, mode="thick"):
    '''Returns an estimte of the accretion fraction facc for a thick or thin
    boundary layer.  Facc is calculated as z/4r, where z is the scale height
    of disk.  If mode = "thick", the disk scale height is assumed to be equal 
    to the disk scale height. If mode = "thin", use approximation of Tylenda
    (1981).  See Patterson & Raymond, 1985 for details.'''
    if mode.lower() == "thick":
        z = (np.sqrt(radius(mwd)) ** 2.0 * \
            (6.96e-4 * (mwd/0.7) ** -0.85 * (mdot * MDOTGS / 1e18) ** -0.22 + \
            (7.29e-4 * (mwd/0.7) ** 0.8 * (mdot * MDOTGS / 1e18))))
        facc = z / 4 / radius(mwd)
        return facc
    elif mode.lower() == "thin":
        z = 3.85e8 * (mwd/0.7) ** 0.1 * (mdot * MDOTGS / 1e15) ** -0.5
        facc = z / 4 / radius(mwd)
        return facc
    else:
        print("mode must be thick or thin")

def coolratio(nebl, urad, tshock):
    ratio = 7.5e-5 * nebl / (urad * tshock**0.5)
    return ratio

mwd = 1.3
mdot = 3e-9
mdotthick = 1e-9
frad = 1.0
mode = "thin"
vacc = 2000

n = nebl(mdot, mwd, facc(mwd, mdot, mode=mode), vacc) * np.exp(-5)
ts = tshock(mwd)
lrad = G * mwd * MSOL * mdotthick * MDOTGS / radius(mwd)
lrad = 4 * np.pi * radius(mwd)**2 * 5.67e-5 * 120000**4.0
print(lrad)
ur = urad(lrad, mwd, frad)
ratio = coolratio(n, ur, ts)
print(ratio)

