from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_library
from JWST_MG.delta_c import delta_c
from JWST_MG.HMF import HMF
from JWST_MG.SMF import SMF

class NaNException(Exception):
    """Integrator hit a NaN."""

    pass

class SMD:
   ########################################################################
    # Initialize a class
    # float a - scale factor value (related to redshift via a = 1/(1+z))
    # string model - model of MG for the derivation of mu parameter
    # string model_H - model of MG for H(a)
    # string model_SFR - model that determines SMHR
    # float par1, par2 - corresponding MG parameters
    # Masses - array of CDM halo masses
    # f0 - parameter for "double-power" SMHR, otherwise None
    ########################################################################

    def __init__(self, a, model, model_H, model_SFR, par1, par2, Masses, f0=None):
        self.a = a
        self.model = model
        self.model_H = model_H
        self.model_SFR = model_SFR
        self.par1 = par1
        self.par2 = par2
        self.Masses = Masses
        self.f0 = f0


    # Taken from https://hmf.readthedocs.io/en/latest/_modules/hmf/mass_function/integrate_hmf.html
    def hmf_integral_gtm(self, M, dndm, mass_density=False):
    
        # Eliminate NaN's
        m = M[np.logical_not(np.isnan(dndm))]
        dndm = dndm[np.logical_not(np.isnan(dndm))]
        dndlnm = m * dndm

        if len(m) < 4:
            raise JWST_MG.SMD.NaNException(
                "There are too few real numbers in dndm: len(dndm) = %s, #NaN's = %s"
                % (len(M), len(M) - len(dndm))
            )

        # Calculate the mass function (and its integral) from the highest M up to 10**18
        if m[-1] < m[0] * 10**18 / m[3]:
            m_upper = np.arange(
                np.log(m[-1]), np.log(10**18), np.log(m[1]) - np.log(m[0])
            )
            mf_func = spline(np.log(m), np.log(dndlnm), k=1)
            mf = mf_func(m_upper)

            if not mass_density:
                int_upper = intg.simps(np.exp(mf), dx=m_upper[2] - m_upper[1], even="first")
            else:
                int_upper = intg.simps(
                    np.exp(m_upper + mf), dx=m_upper[2] - m_upper[1], even="first"
                )
        else:
            int_upper = 0

        # Calculate the cumulative integral (backwards) of [m*]dndlnm
        if not mass_density:
            ngtm = np.concatenate(
                (
                    intg.cumtrapz(dndlnm[::-1], dx=np.log(m[1]) - np.log(m[0]))[::-1],
                    np.zeros(1),
                )
            )
        else:
            ngtm = np.concatenate(
                (
                    intg.cumtrapz(m[::-1] * dndlnm[::-1], dx=np.log(m[1]) - np.log(m[0]))[
                        ::-1
                    ],
                    np.zeros(1),
                )
            )

        return ngtm + int_upper

    def SMD(self, Masses, rhoM, a, model_H, model,model_SFR, par1, par2, f0):
        HMF_library = HMF(a, model, model_H, par1, par2, Masses/h)
        HMF_fid = HMF_library.ST_mass_function(rhoM, Masses/h, a, model_H, model, par1, par2)
        SMF_library = SMF(a, model, model_H, model_SFR, par1, par2, Masses, f0)
        fstar = SMF_library.epsilon(Masses, model_SFR, a, f0)*Omegab0/Omegam0
        print(fstar)
        SMD_fid = fstar*self.hmf_integral_gtm(Masses, HMF_fid, mass_density = True)
        return fstar*Masses, SMD_fid