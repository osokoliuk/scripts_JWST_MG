import numpy as np
import os

def SMF_z4_obs():
    x = []
    y = []
    yerr_up = []
    yerr_down = []
    path = os.path.dirname(os.path.realpath(__file__))

    Duncan = np.loadtxt(path + "/Duncan_z4.dat")
    np.concatenate((x, 10**Duncan[:,0]), axis=None)
    np.concatenate((x, Duncan[:,1]), axis=None)
    np.concatenate((yerr_up, Duncan[:,2]), axis=None)
    np.concatenate((yerr_down, Duncan[:,3]), axis=None)

    Song = np.loadtxt(path + "/Song_z4.dat")
    np.concatenate((x, 10**Song[:,0]), axis=None)
    np.concatenate((x, 10**Song[:,1]), axis=None)
    np.concatenate((yerr_up, 10**Song[:,2]), axis=None)
    np.concatenate((yerr_down, 10**Song[:,3]), axis=None)

    Navarro = np.loadtxt(path + "/Navarro_z4.dat")
    np.concatenate((x, 10**Navarro[:,0]), axis=None)
    np.concatenate((x, 1e-4*Navarro[:,1]), axis=None)
    np.concatenate((yerr_up, 1e-4*Navarro[:,2]), axis=None)
    np.concatenate((yerr_down, 1e-4*Navarro[:,2]), axis=None)

    return x, y, yerr_down, yerr_up


def SMF_z4_obs_dict():
    x = {}
    y = {}
    yerr_up = {}
    yerr_down = {}
    path = os.path.dirname(os.path.realpath(__file__))
    Duncan = np.loadtxt(path + "/Duncan_z4.dat")
    x.update({'Duncan': 10**Duncan[:,0]})
    y.update({'Duncan': Duncan[:,1]})
    yerr_down.update({'Duncan': Duncan[:,2]})
    yerr_up.update({'Duncan': Duncan[:,3]})

    Song = np.loadtxt(path + "/Song_z4.dat")
    x.update({'Song': 10**Song[:,0]})
    y.update({'Song': 10**Song[:,1]})
    yerr_down.update({'Song': 10**Song[:,2]})
    yerr_up.update({'Song': 10**Song[:,3]})
    
    Navarro = np.loadtxt(path + "/Navarro_z4.dat")
    x.update({'Navarro': 10**Navarro[:,0]})
    y.update({'Navarro': 1e-4*Navarro[:,1]})
    yerr_down.update({'Navarro': 1e-4*Navarro[:,2]})
    yerr_up.update({'Navarro': 1e-4*Navarro[:,2]})

    return x, y, yerr_down, yerr_up

