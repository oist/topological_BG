#!/usr/bin/env python 

bgParams =\
{
    "FSI_iaf": {
        "C_m": 3.1,
        "V_th": 16.0,
        "tau_m": 3.1
    },
    "GFSI": 1.0,
    "GGPe": 1.0,
    "GGPe_STN": 1.0,
    "GGPi": 1.0,
    "GMSN": 1.0,
    "GMSN_MSN": 1.0, #3.5,
    "GMSN_GPx": 0.3,  #to slow down the inhibition on GPi and GPe from MSNs  (only used in case of dopamine-based plasticity on CS+ episodes)
    "GPe_iaf": {
        "C_m": 14.0,
        "V_th": 11.0,
        "tau_m": 14.0
    },
    "GPi_iaf": {
        "C_m": 14.0,
        "V_th": 6.0,
        "tau_m": 14.0
    },
    "GSTN": 1.0,
    "GSTN_GPi": 1.0,
    "IeFSI": 8.0,
    "IeGPe": 11.0,
    "IeGPi": 8.5,
    "IeMSN": 24.5, #26.0, #24.5,#26.0,
    "IeSTN": 9.0,
    "MSN_iaf": {
        "C_m": 13.0,
        "V_th": 30.0,
        "tau_m": 13.0
    },
    "ProjPercent": {
        "CMPf->FSI": 1.0,
        "CMPf->GPe": 1.0,
        "CMPf->GPi": 1.0,
        "CMPf->MSN": 1.0,
        "CMPf->STN": 1.0,
        "CSN->FSI": 1.0,
        "CSN->MSN": 1.0,
        "FSI->FSI": 1.0,
        "FSI->MSN": 1.0,
        "GPe->FSI": 0.16,
        "GPe->GPe": 0.84,
        "GPe->GPi": 0.84,
        "GPe->MSN": 0.16,
        "GPe->STN": 1.0,
        "MSN->GPe": 1.0,
        "MSN->GPi": 0.82,
        "MSN->MSN": 1.0,
        "PTN->FSI": 1.0,
        "PTN->MSN": 1.0,
        "PTN->STN": 1.0,
        "STN->FSI": 0.17,
        "STN->GPe": 0.83,
        "STN->GPi": 0.72,
        "STN->MSN": 0.17
    },
    "RedundancyType": "outDegreeAbs",
    "Ri": 2.0,
    "Rm": 2.0,
    "STN_iaf": {
        "C_m": 6.0,
        "V_th": 26.0,
        "tau_m": 6.0
    },
    "alpha": {
        "CMPf->FSI": 1053,
        "CMPf->GPe": 79,
        "CMPf->GPi": 131,
        "CMPf->MSN": 4965,
        "CMPf->STN": 76,
        "CSN->FSI": 250,
        "CSN->MSN": 342,
        "FSI->FSI": 116,
        "FSI->MSN": 4362,
        "GPe->FSI": 353,
        "GPe->GPe": 38,
        "GPe->GPi": 16,
        "GPe->MSN": 0,
        "GPe->STN": 19,
        "MSN->GPe": 171,
        "MSN->GPi": 210,
        "MSN->MSN": 210,
        "PTN->FSI": 5,
        "PTN->MSN": 5,
        "PTN->STN": 259,
        "STN->FSI": 91,
        "STN->GPe": 428,
        "STN->GPi": 233,
        "STN->MSN": 0
    },
    "asymmetry_1": 0.5,  #30%, structural asymmetry. this parameter is kept fixed. (MSN_D1->MSN_D1, MSN_D2->MSN_D2, MSN_D2->MSN_D1)
    "asymmetry_2": 0.1, #6%, structural asummetry. we explore this asymmetry (MSN_D1->MSN_D2)
    "MSND2_mod":  1.0, # 1=no modulation. modulation of synapses from MSND2->MSNs (used only during dopamine dips, on CS- episodes) 
    "syn_asymm": 2.0, # functional asymmetry. This will increase the synaptic streght from MSND2 to others MSNs.
    "cTypeCMPfFSI": "diffuse",
    "cTypeCMPfGPe": "diffuse",
    "cTypeCMPfGPi": "diffuse",
    "cTypeCMPfMSN": "diffuse",
    "cTypeCMPfSTN": "diffuse",
    "cTypeCSNFSI": "focused",
    "cTypeCSNMSN": "focused",
    "cTypeFSIFSI": "diffuse",
    "cTypeFSIMSN": "diffuse",
    "cTypeGPeFSI": "diffuse",
    "cTypeGPeGPe": "diffuse",
    "cTypeGPeGPi": "focused",
    "cTypeGPeMSN": "diffuse",
    "cTypeGPeSTN": "focused",
    "cTypeMSNGPe": "focused",
    "cTypeMSNGPi": "focused",
    "cTypeMSNMSN": "focused",
    "cTypePTNFSI": "focused",
    "cTypePTNMSN": "focused",
    "cTypePTNSTN": "focused",
    "cTypeSTNFSI": "diffuse",
    "cTypeSTNGPe": "diffuse",
    "cTypeSTNGPi": "diffuse",
    "cTypeSTNMSN": "diffuse",
    "channels": True,
    "circle_center": [],
    "common_iaf": {
        "E_L": 0.0,
        "I_e": 0.0,
        "V_m": 0.0,
        "V_min": -20.0,
        "V_reset": 0.0,
        "V_th": 10.0,
        "t_ref": 2.0,
        "tau_syn": [
            1.8393972058572117,
            36.787944117144235,
            1.8393972058572117
        ]
    },
    "countCMPf": 86000.0,
    "countCSN": None,
    "countFSI": 532000.0,
    "countGPe": 251000.0,
    "countGPi": 143000.0,
    "countMSN": 26448000.0,
    "countPTN": None,
    "countSTN": 77000.0,
    "distcontact": {
        "CMPf->FSI": 0.06,
        "CMPf->GPe": 0.0,
        "CMPf->GPi": 0.48,
        "CMPf->MSN": 0.27,
        "CMPf->STN": 0.46,
        "CSN->FSI": 0.82,
        "CSN->MSN": 0.95,
        "FSI->FSI": 0.16,
        "FSI->MSN": 0.19,
        "GPe->FSI": 0.58,
        "GPe->GPe": 0.01,
        "GPe->GPi": 0.13,
        "GPe->MSN": 0.06,
        "GPe->STN": 0.58,
        "MSN->GPe": 0.48,
        "MSN->GPi": 0.59,
        "MSN->MSN": 0.77,
        "PTN->FSI": 0.7,
        "PTN->MSN": 0.98,
        "PTN->STN": 0.97,
        "STN->FSI": 0.41,
        "STN->GPe": 0.3,
        "STN->GPi": 0.59,
        "STN->MSN": 0.16
    },
    "dx": {
        "FSI": 1.5e-06,
        "GPe": 1.7e-06,
        "GPi": 1.2e-06,
        "MSN": 1e-06,
        "STN": 1.5e-06
    },
    "lx": {
        "FSI": 0.000961,
        "GPe": 0.000865,
        "GPi": 0.001132,
        "MSN": 0.000619,
        "STN": 0.00075
    },
    "nbCMPf": 36000.0,
    "nbCSN": 36000.0,
    "nbFSI": 636.0,
    "nbGPe": 300.0,
    "nbGPi": 168.0,
    "nbMSN": 31728.0,
    "nbPTN": 36000.0,
    "nbSTN": 96.0,
    "normalrate": {
        "CMPf": [
            4.0,
            34.0
        ],
        "CSN": [
            2.0,
            19.7
        ],
        "FSI": [
            7.8,
            14.0
        ],
        "GPe": [
            55.7,
            74.5
        ],
        "GPi": [
            59.1,
            79.5
        ],
        "MSN": [
            0.05,
            1
        ],
        "PTN": [
            15.0,
            46.3
        ],
        "STN": [
            15.2,
            22.8
        ]
    },
    "num_neurons": 1000,
    "overlap_d1d2": 0.1, #direct and indirect pathways overlap.
    "parrotCMPf": True,
    "plast_gain": 0.65,
    "plastic_syn": True, #needed for stdp for cortical afferents into the striatum
    "redundancyCMPfFSI": 3,
    "redundancyCMPfGPe": 3,
    "redundancyCMPfGPi": 3,
    "redundancyCMPfMSN": 3,
    "redundancyCMPfSTN": 3,
    "redundancyCSNFSI": 3,
    "redundancyCSNMSN": 3,
    "redundancyFSIFSI": 3,
    "redundancyFSIMSN": 3,
    "redundancyGPeFSI": 3,
    "redundancyGPeGPe": 3,
    "redundancyGPeGPi": 3,
    "redundancyGPeMSN": 3,
    "redundancyGPeSTN": 3,
    "redundancyMSNGPe": 3,
    "redundancyMSNGPi": 3,
    "redundancyMSNMSN": 3,
    "redundancyPTNFSI": 3,
    "redundancyPTNMSN": 3,
    "redundancyPTNSTN": 3,
    "redundancySTNFSI": 3,
    "redundancySTNGPe": 3,
    "redundancySTNGPi": 3,
    "redundancySTNMSN": 3,
    "spread_diffuse": 2.0,
    "spread_focused": 0.15,
    "stochastic_delays": None,
    "tau": {
        "CMPf->FSI": 7.0,
        "CMPf->GPe": 7.0,
        "CMPf->GPi": 7.0,
        "CMPf->MSN": 7.0,
        "CMPf->STN": 7.0,
        "CSN->FSI": 7.0,
        "CSN->MSN": 7.0,
        "FSI->FSI": 1.0,
        "FSI->MSN": 1.0,
        "GPe->FSI": 3.0,
        "GPe->GPe": 1.0,
        "GPe->GPi": 3.0,
        "GPe->MSN": 3.0,
        "GPe->STN": 10.0,
        "MSN->GPe": 7.0,
        "MSN->GPi": 11.0,
        "MSN->MSN": 1.0,
        "PTN->FSI": 3.0,
        "PTN->MSN": 3.0,
        "PTN->STN": 3.0,
        "STN->FSI": 3.0,
        "STN->GPe": 3.0,
        "STN->GPi": 3.0,
        "STN->MSN": 3.0
    },
    "wPSP": [
        1.0,
        0.025,
        -0.25
    ]
}
