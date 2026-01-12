#!/usr/bin/env python
# -*- coding: utf-8 -*-

##
## fetch_params.py
##
## This file contains routines that set the parameter files and return an object containing everything required by the instantiate_XYZ functions of ini_all.py
# ideally, functions below should read from parameter files, but it may be easier to hardcode values here as a first try

def read_sim():
  try:
    from runpy import run_path
    file_params = run_path('simParams.py', init_globals=globals())
    sim_params = file_params['simParams']
    return sim_params
  except:
    raise ImportError('The simulation parameters could not be loaded. Please make sure that the file `simParams.py` exists and is a valid python defining the variable "simParams".')

def read_bg():
  try:
    from runpy import run_path
    file_params = run_path('bgParams.py', init_globals=globals())
    bg_params = file_params['bgParams']
    return bg_params
  except:
    raise ImportError('The BG-region parameters could not be loaded. Please make sure that the file `bgParams.py` exists and is a valid python defining the variable "bgParams".')


