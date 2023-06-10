%------------------------------------------------------------------------------
%
%                          LAMBERBATTIN
%
%   this subroutine solves Lambert's problem using Battins method. The method
%   is developed in Battin (1987). It uses contiNued fractions to speed the
%   solution and has several parameters that are defined differently than
%   the traditional Gaussian technique.
%
% Inputs:         Description                    Range/Units
%   ro          - IJK Position vector 1          m
%   r           - IJK Position vector 2          m
%   dm          - direction of motion            'pro','retro'
%   Dtsec       - Time between ro and r          s
%
% OutPuts:
%   vo          - IJK Velocity vector            m/s
%   v           - IJK Velocity vector            m/s
%
% Reference:
% Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill
% , New York; 3rd edition(2007).
% 
% Last modified:   2015/08/12   M. Mahooti
% 
%------------------------------------------------------------------------------ 
import numpy as np

def LAMBERTBATTIN(ro, r, dm, Dtsec):
    
    small = 0.000001
    mu = 3.986004418e14   # m3/s2
    y1 = 0
    magr = np.linalg.norm(r)
    magro = np.linalg.norm(ro)
    CosDeltaNu = np.dot(ro, r)/(magro*magr)
    rcrossr = np.cross(ro, r)
    magrcrossr = np.linalg.norm(rcrossr)
    if dm == 'pro':
        SinDeltaNu = magrcrossr/(magro*magr)
    else:
        SinDeltaNu = -magrcrossr/(magro*magr)
    DNu = np.arctan2(SinDeltaNu, CosDeltaNu)

    # the angle needs to be positive to work for the long way
    if DNu < 0.0:
        DNu = 2.0*np.pi + DNu

    RoR   = magr/magro
    eps   = RoR - 1.0
    tan2w = 0.25*eps*eps / (np.sqrt(RoR) + RoR*(2.0 + np.sqrt(RoR)))
    rp    = np.sqrt(magro*magr)*((np.cos(DNu*0.25))**2 + tan2w)
    if (DNu < np.pi):
        L = ((np.sin(DNu*0.25))**2 + tan2w) / ((np.sin(DNu*0.25))**2 + tan2w + np.cos(DNu*0.5))
    else:
        L = ((np.cos(DNu*0.25))**2 + tan2w - np.cos(DNu*0.5)) / ((np.cos(DNu*0.25))**2 + tan2w)
    m    = mu*Dtsec*Dtsec / (8.0*rp*rp*rp)
    x    = 10.0
    xn   = L
    chord = np.sqrt(magro*magro + magr*magr - 2.0*magro*magr*np.cos(DNu))
    s    = (magro + magr + chord)*0.5
    lim1 = np.sqrt(m/L)
    
    Loops = 1
    while True:
    x = xn
    tempx = seebatt(x)
    Denom = 1.0 / ((1.0+2.0*x+L) * (4.0*x + tempx*(3.0+x)))
    h1 = (L+x)**2 * (1.0+3.0*x+tempx) * Denom
    h2 = m*(x-L+tempx) * Denom

    # ----------------------- Evaluate CUBIC ------------------
    b = 0.25*27.0*h2 / ((1.0+h1)**3)
    if b < -1.0: # reset the initial condition
        xn = 1.0 - 2.0*l
    else:
        if y1 > lim1:
            xn = xn * (lim1/y1)
        else:
            u = 0.5*b / (1.0 + sqrt(1.0 + b))
            k2 = seebattk(u)
            y = ((1.0+h1)/3.0) * (2.0 + sqrt(1.0+b)/(1.0+2.0*u*k2*k2))
            xn = sqrt(((1.0-L)*0.5)**2 + m/(y*y)) - (1.0+L)*0.5
    Loops = Loops + 1
    y1 = sqrt(m/((L+x)*(1.0+x)))
    if abs(xn-x) < small and Loops > 30:
        break

    a = mu*Dtsec*Dtsec / (16.0*rp*rp*xn*y*y)
	
    # ------------------ Find Eccentric anomalies -----------------
# ------------------------ Hyperbolic -------------------------
if a < -small:
    arg1 = np.sqrt(s / (-2.0 * a))
    arg2 = np.sqrt((s - chord) / (-2.0 * a))
    # ------- Evaluate f and g functions --------
    AlpH = 2.0 * np.arcsinh(arg1)
    BetH = 2.0 * np.arcsinh(arg2)
    DH = AlpH - BetH
    F = 1.0 - (a / magro) * (1.0 - np.cosh(DH))
    GDot = 1.0 - (a / magr) * (1.0 - np.cosh(DH))
    G = Dtsec - np.sqrt(-a * a * a / mu) * (np.sinh(DH) - DH)
else:
    # ------------------------ Elliptical ---------------------
    if a > small:
        arg1 = np.sqrt(s / (2.0 * a))
        arg2 = np.sqrt((s - chord) / (2.0 * a))
        Sinv = arg2
        Cosv = np.sqrt(1.0 - (magro + magr - chord) / (4.0 * a))
        BetE = 2.0 * np.arccos(Cosv)
        BetE = 2.0 * np.arcsin(Sinv)
        if DNu > np.pi:
            BetE = -BetE
        Cosv = np.sqrt(1.0 - s / (2.0 * a))
        Sinv = arg1
        am = s * 0.5
        ae = np.pi
        be = 2.0 * np.arcsin(np.sqrt((s - chord) / s))
        tm = np.sqrt(am * am * am / mu) * (ae - (be - np.sin(be)))
        if Dtsec > tm:
            AlpE = 2.0 * np.pi - 2.0 * np.arcsin(Sinv)
        else:
            AlpE = 2.0 * np.arcsin(Sinv)
        DE = AlpE - BetE
        F = 1.0 - (a / magro) * (1.0 - np.cos(DE))
        GDot = 1.0 - (a / magr) * (1.0 - np.cos(DE))
        G = Dtsec - np.sqrt(a * a * a / mu) * (DE - np.sin(DE))
    else:
        # --------------------- Parabolic ---------------------
        arg1 = 0.0
        arg2 = 0.0
        raise ValueError('a parabolic orbit')

    for i in range(3):
    vo[i] = (r[i] - F*ro[i]) / G
    v[i] = (GDot*r[i] - ro[i]) / G

    return vo, v