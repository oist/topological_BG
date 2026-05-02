#!/usr/bin/env python
# -*- coding: utf-8 -*-

##
## ini_all.py
##
## This contains region-specific initialization functions, based on provided parameters.

# nest_routines contains all nest-specific calls
import numpy as np
import nest_routine
import nest
import time
import collections

# the functions below contain the logics of the instantiation, without any direct call to nest



def instantiate_bg(bg_params, fake_inputs=False, ctx_inputs=None,scalefactor=[1,1]):
  bg_layers = {}
  ctx_bg_input = {}
  #update bg_params with the scaled number of neurons
  # 1 mmmm2 is the base for # of neurons, if the surface increases more neurons are added
  for nucleus in ['MSN','FSI','STN','GPe','GPi','CSN','PTN','CMPf']:
    bg_params['nb' + nucleus] = bg_params['nb' + nucleus] * scalefactor[0] * scalefactor[1]
  
  print('###########################################')
  print("BG instantiation")
  for nucleus in ['GPi','MSN','FSI','STN','GPe','GPi_fake']:
    print("Creating MAIN nucleus "+nucleus+"...")
    result = nest_routine.create_layers_bg(bg_params, nucleus, scalefactor=scalefactor)
    if nucleus == 'MSN':
      bg_layers['MSN_d1'], bg_layers['MSN_d2'] = result
    else:
      bg_layers[nucleus] = result
  
  #connect GPi2D to GPi3D(fake)
  print('connecting GPi 2D to fake GPi 3D ...')
  nest_routine.connect_GPi2d_GPi3d(bg_layers['GPi'], bg_layers['GPi_fake'])
  
  print('creating fake input layers in the BG ... ')
  rate=0.
  if fake_inputs:
    for fake_nucleus in ['CSN','PTN','CMPf']:
      rate = bg_params['normalrate'][fake_nucleus][0]      
      mirror_neurons = None
      print("Creating fake input "+fake_nucleus+" with fire rate "+str(rate)+" Hz...")
      bg_layers[fake_nucleus]= nest_routine.create_layers_bg(bg_params, fake_nucleus, fake=rate, mirror_neurons=mirror_neurons,mirror_pos=None,scalefactor=scalefactor)

  ###################################################################################################
  ##### creation volume transmitter and initial setting for STD-Dopa synapse #########################
  ##### this synapses will behave as static synapses unless dopamine is released (reward)    #########
  ###################################################################################################
  if bg_params['plastic_syn']:
    bg_params['vt_d1']=nest.Create('volume_transmitter')
    bg_params['vt_d2']=nest.Create('volume_transmitter')

    # NEST 3: parameter renamed from 'vt' to 'volume_transmitter'; NodeCollection, set via SetDefaults
    nest.CopyModel('stdp_dopamine_synapse_lbl', 'syn_d1')
    bg_params['wr_d1'] = nest.Create('weight_recorder', 1, {'record_to': 'ascii'})
    nest.SetDefaults('syn_d1', {'volume_transmitter': bg_params['vt_d1'],
                                'A_plus': 0.013, 'A_minus': (0.013/4.0), 'Wmax': 4., 'b': 0., 'n': 0., 'c': 0.,
                                'tau_plus': 20., 'tau_n': 100.0, 'tau_c': 700.})

    nest.CopyModel('stdp_dopamine_synapse_lbl', 'syn_d2')
    bg_params['wr_d2'] = nest.Create('weight_recorder', 1, {'record_to': 'ascii'})
    nest.SetDefaults('syn_d2', {'volume_transmitter': bg_params['vt_d2'],
                                'A_plus': 0.013, 'A_minus': -0.013, 'Wmax': 4., 'b': 0., 'n': 0., 'c': 0.,
                                'tau_plus': 20., 'tau_n': 100.0, 'tau_c': 700.})


    with open('./log/weights_d1.txt','a') as file:
      file.write('initial setting for stdp_dopamine_synapse_lbl d1: '+str(nest.GetDefaults('syn_d1'))+' \n')
    with open('./log/weights_d2.txt','a') as file:
      file.write('initial setting for stdp_dopamine_synapse_lbl d2: '+str(nest.GetDefaults('syn_d2'))+' \n')


  print("BG connect layers")
  #########################
  bg_params['alpha'] = collections.OrderedDict(sorted(bg_params['alpha'].items(), key=lambda t: t[0]))
  #########################
  for connection in bg_params['alpha']:
    src = connection[0:3]
    tgt = connection[-3:]
    if src == 'CMP':
      src = 'CMPf' # string split on '->' would be better
    # if not fake_inputs, wire only in case of intra-BG connection
    if fake_inputs or src in ['MSN','FSI','STN','GPe','GPi']:
      if src in ['MSN','FSI','GPe','GPi']:
        nType = 'in'
      else:
        nType = 'ex'
      nest_routine.connect_layers_bg(bg_params, nType, bg_layers, src, tgt, projType=bg_params['cType'+src+tgt], redundancy=bg_params['redundancy'+src+tgt], RedundancyType=bg_params['RedundancyType'],verbose=True,scalefactor=scalefactor)
    #### UNCOMMENT LINES BELOW IF FILES FOR BG CONNECTIONS ARE NEEDED #####
    #if tgt=='MSN' or src=='MSN': #added to save connections to MSN as txt files
    #  nest_routine.get_connections(src,tgt,bg_layers[src],bg_layers[tgt])
  return bg_layers,ctx_bg_input    


