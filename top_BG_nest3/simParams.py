#!/usr/bin/env python 
# this file set up parameters for the simulation.


import random

simParams =\
{
    "channels": True,
    "channels_nb": 6,
    "channels_radius": 0.12,
    "circle_center": [],
    "dt": "0.1",
    "hex_radius": 0.24,
    "initial_ignore": 0.0,
    "msd": random.randint(0, 1000), 
    "nbcpu": 10,
    "nbnodes": 1,
    "overwrite_files": True,
    "scalefactor": [
        1.0,
        1.0
    ],
    "simDuration": 3000.0,
    "sim_model": {
        "action_selection": {
            "on": False,
            "regions": {
                "BG": True,
                "CB_M1": False,
                "CB_S1": False,
                "M1": False,
                "S1": False,
                "TH_M1": False,
                "TH_S1": False
            }
        },
        "plasticity": {
            "on": True,
            "regions": {
                "BG": True,
                "CB_M1": False,
                "CB_S1": False,
                "M1": False,
                "S1": False,
                "TH_M1": False,
                "TH_S1": False
            }
        },
        "resting_state": {
            "on": False,
            "regions": {
                "BG": True,
                "CB_M1": False,
                "CB_S1": False,
                "M1": False,
                "S1": False,
                "TH_M1": False,
                "TH_S1": False
            }
        },
        
    },
    "start_time_sp": 1000.0,
    "whichSim": "stim_all_model"
}
