#!/usr/bin/env python
# -*- coding: utf-8 -*-

##
## nest_routine.py
##
## This script defines creation, connection, and simulation routines using PyNest
##
##

import nest.topology as ntop
from nest.lib.hl_api_info import SetStatus
import nest

import numpy as np
import math
import time
import collections
import os
#import pandas as pd
import random
import itertools
###########
# General #
###########
#### global variable to be updated when initializing nest #######
pyrngs = []


#-------------------------------------------------------------------------------
# Nest initialization
#-------------------------------------------------------------------------------
def initialize_nest(sim_params):
  nest.set_verbosity("M_WARNING")
  nest.SetKernelStatus({"overwrite_files": sim_params['overwrite_files']}) # should we erase previous traces when redoing a simulation?
  nest.SetKernelStatus({'local_num_threads': int(sim_params['nbcpu'])})
  nest.SetKernelStatus({"data_path": 'log'})
  if sim_params['dt'] != '0.1':
    nest.SetKernelStatus({'resolution': float(sim_params['dt'])})
  ####### adding changes to nest seeds for independent experiments ##########
  ### changing python seeds ####
  N_vp = nest.GetKernelStatus(['total_num_virtual_procs'])[0]
  print('Number of virtual processes: ',N_vp)
  global pyrngs
  pyrngs = [np.random.RandomState(s) for s in range(sim_params['msd'], sim_params['msd']+N_vp)]
  ### global nest rng ###
  nest.SetKernelStatus({'grng_seed' : sim_params['msd']+N_vp})
  ### per process rng #####
  nest.SetKernelStatus({'rng_seeds' : range(sim_params['msd']+N_vp+1, sim_params['msd']+2*N_vp+1)})


#-------------------------------------------------------------------------------
# Instantiate a spike detector and connects it to the entire layer `layer_gid`
#-------------------------------------------------------------------------------
def layer_spike_detector(layer_gid, layer_name, ignore_time, params={"withgid": True, "withtime": True, "to_file": True, "fbuffer_size":8192}):
#def layer_spike_detector(layer_gid, layer_name, params={"withgid": True, "withtime": True, "to_file": True, 'fbuffer_size': 8192}):
    print ('spike detector for '+layer_name)
    params.update({'label': layer_name, "start": float(ignore_time)})
    #params.update({'label': layer_name})
    detector = nest.Create("spike_detector", params= params)
    nest.Connect(pre=nest.GetNodes(layer_gid)[0], post=detector)
    return detector

def average_fr(detector, simDuration, n):
    return nest.GetStatus(detector, 'n_events')[0] / (float(simDuration) * float(n) / 1000.)

def instantaneous_fr(detector,sim_start,sim_end,n,dt=0.1):
    At = []
    spikes = nest.GetStatus(detector)[0]['events']
    for i in np.arange(int(sim_start),int(sim_end),dt): #np.arange(int(simtime))
      a = ((i <= spikes['times']) & (spikes['times'] < (i+1))).sum()        
      At.append(a)
    At = np.array(At)/float(n)*1000.  # empirical instantaneous activity At
    return At


#-------------------------------------------------------------------------------
# Returns the number of neurons inside a layer
#-------------------------------------------------------------------------------
def count_layer(layer_gid):
  return len(nest.GetNodes(layer_gid)[0])

#-------------------------------------------------------------------------------
# Returns the connections of neurons inside a layer -sun-20180912
#-------------------------------------------------------------------------------
def get_connection(gids):
  return nest.GetConnections(gids)


def get_columns_data(layer_name,circle_center,radius_small,my_area=None):
    #### example ##########
    #pkj_M1 = get_columns_data('CB_M1_layer_pkj',bg_params['circle_center'] ,bg_params['channels_radius'])
    ##########################
    gid_pos = np.loadtxt('./log/'+layer_name+'.txt')
    print('gids and pos  ', len(gid_pos))
    circle_gids = []
    for i in np.arange(len(circle_center)):
        idx = np.where(np.linalg.norm([(gid_pos[:,1]-circle_center[i][0]),(gid_pos[:,2]-circle_center[i][1])],axis=0)<=radius_small)[0]
        print('number of neurons in channel ',str(i),'for ',layer_name,' : ',str(len(idx)))
        circle_gids.append([[int(x[0]),x[1:].tolist()] for x in gid_pos[idx,:]])
    return circle_gids #circle_gids


AMPASynapseCounter_bg = 0 # initialize global counter variable for AMPA/NMDA colocalization in BG (unfortunate design choice, but required by nest fast connect procedure)

## Function hex_corner provides vertex coordinates of a hexagon, given a center, a radius (size) and the vertex id.
def hex_corner(center,size,i):
    angle_deg = 60 * i - 30
    angle_rad = np.pi / 180 * angle_deg
    return [center[0] + size * np.cos(angle_rad),center[1] + size * np.sin(angle_rad)]


#define the centers that will connect ctx to bg, and store them at bg_params['circle_center']
#centers must be within grid 2D dimensions.
def get_channel_centers(bg_params,hex_center=[0,0],ci=6,hex_radius=0.240):
    center_aux = []
    if bg_params['channels']:
        if len(bg_params['circle_center'])==0: #must be done before bg instantiation.
            for i in np.arange(ci):
                x_y = hex_corner(hex_center,hex_radius,i) #center, radius, vertex id # gives x,y of an hexagon vertexs.
                center_aux.append(x_y)
                            #bg_params['circle_center'].append(x_y)
            np.savetxt('./log/centers.txt',center_aux) #save the centers.
            print('generated centers: ',center_aux)
    return center_aux


def apply_fixed_stimulus(time_start, time_stop, psg, stimulus):
    nest.SetStatus([psg],{'rate':stimulus,'start':time_start,'stop':time_stop})


def get_targets_mean_rates(columns_gids,spike_detector,idx_dt,t_start,t_end):
    targets_rates = []
    #get the 3 top GIDs:
    for col in np.arange(5):#[1,2,3]: #indexes of the top 3 columns
        gid_l = [gid[0] for gid in columns_gids[col]]
        mask = np.in1d(spike_detector['senders'][idx_dt],gid_l)
        mean_act = len(spike_detector['times'][idx_dt][mask])/((t_end-t_start)*float(len(columns_gids[col])))*1000. 
        #mean_act = len(spike_detector['times'][idx_dt][mask])/((dt_task-2.*dt_task_no_stim)*float(len(columns_gids[col])))*1000. 
        print('mean firing rate for col ',col,':  ',mean_act)
        targets_rates.append(mean_act)
    return targets_rates
    

def create_psg_channels(syn_weight, syn_delay, channels_nb):
    psg = nest.Create('poisson_generator', channels_nb)
    syn = {'weight':syn_weight, 'delay': syn_delay}  # syn_PG_S1_L5A_Pyr={'weight': 2.0,'delay': 1.5}
    return psg, syn


def connect_psg_to_channels(columns_gids, psg, syn):
    circle_j_gids_nb = {}
    for j in np.arange(len(columns_gids)):  # for each circle of gids
        circle_j_gids = [k[0] for k in columns_gids[j]]
        print('circle  ', j, ' contains #: ', len(circle_j_gids))
        circle_j_gids_nb[str(j)] = len(circle_j_gids)
        nest.Connect(pre=[psg[j]], post=circle_j_gids, syn_spec=syn)



### function to define BG grid positions in 2D
### parameters: 
# nbCh: number of channels (always 1)
# sim_pts: number of points to generate
# a0, a1:  -x shift, distance from starting point (x axis)
# b0, b1:  -y shift, distance from starting point (y axis)
# -----------------------------------------------------
def grid_positions(nbCh, sim_pts,a0,a1,b0,b1):
    #N_vp = nest.GetKernelStatus(['total_num_virtual_procs'])[0]
    #pyrngs = [np.random.RandomState(123456)]# for s in ([123456]*N_vp)]
    n = int(sim_pts*nbCh)
    n_squared = np.ceil(np.sqrt(n))
    coord = [[x/n_squared*a1-a0, y/n_squared*b1-b0] for x in np.arange(0,n_squared, dtype=float) for y in np.arange(0,n_squared, dtype=float)]
    # too many points due to square root rounding? remove at random # same random numbers over multiple nodes
    if len(coord) > n:
        coord = np.array(coord)[np.sort(pyrngs[0].choice(range(len(coord)), size=n, replace=False))].tolist()
    aux_x = [coord[i][0] for i in range(len(coord))]
    aux_y = [coord[i][1] for i in range(len(coord))]
    return [aux_x, aux_y]


def grid_uniform_positions(number_of_neurons):
    pos = [[np.random.uniform(-0.5, 0.5), np.random.uniform(-0.5, 0.5)] for j in range(number_of_neurons)]
    return pos

def save_layers_position(layer_name, layer_gid, positions):
    gid_and_positions=np.column_stack((np.array(nest.GetNodes(layer_gid)[0]),positions))
    #if not os.path.exists('log/'+layer_name+'.txt'):
    np.savetxt('log/'+layer_name+'.txt', gid_and_positions, fmt='%1.3f')

    if layer_name=='MSN':
      np.savetxt('log/'+layer_name+'_d1.txt', gid_and_positions[:int(len(gid_and_positions)/2)], fmt='%1.3f')
      np.savetxt('log/'+layer_name+'_d2.txt', gid_and_positions[int(len(gid_and_positions)/2):], fmt='%1.3f')


#-------------------------------------------------------------------------------
# Establishes a topological layer and returns it
# bg_params: basal ganglia parameters
# nucleus: name of the nucleus to instantiate
# fake: numerical value - if 0, then a real population of iaf is instantiated
#                       - if fake > 0, then a Poisson generator population firing at `fake` Hz is instantiated
# force_pop_size: if defined, initialize only this number of neurons
#                 -> this is useful for the cortical connections, as some inputs will be derived from L5A and L5B layers
#-------------------------------------------------------------------------------
def create_layers_bg(bg_params, nucleus, fake=0, mirror_neurons=None, mirror_pos=None,scalefactor=[1,1]):
  
  #define extent and center for 2D layer
  my_extent = [1.*int(scalefactor[0])+1.,1.*int(scalefactor[1])+1.]
  my_center = [0.0, 0.0]

  if mirror_neurons is None:
    # normal case: full input layer is created
    if nucleus =='GPi_fake':
      nucleus_tmp = nucleus[:3]
      pop_size = int(bg_params['nb'+nucleus_tmp])
    else:
      pop_size = int(bg_params['nb' + nucleus]) 
  else:
    # inputs come from existing ctx layer: only a fraction of poisson generators are created
    pop_size = int(bg_params['nb' + nucleus]) - len(mirror_neurons)
  
  print('population size for '+nucleus+': '+str(pop_size))

  if nucleus=='GPi_fake':
    # get xy position from real GPi and add z value
    positions_z = pyrngs[0].uniform(0., 0.5, pop_size).tolist()
    positions = np.loadtxt('./log/'+nucleus[:3]+'.txt') # retrive positions x,y from GPi
    position_nD = [[positions[i][1], positions[i][2], positions_z[i]] for i in range(len(positions))]
    #define extent and center for 3D layer
    my_extent = my_extent + [1.]
    my_center = my_center + [0.]
    print('positions GPi_fake: ',position_nD[:10])
  else:
    if nucleus=='MSN':
      positions = grid_uniform_positions(pop_size) #random positions.
      position_nD = positions
    else:
      if (nucleus=='GPi' or nucleus=='STN'):
        positions = grid_positions(1,pop_size,0.4*scalefactor[0],scalefactor[0]-0.1,0.4*scalefactor[1],scalefactor[1]-0.1)
      else:
        positions = grid_positions(1,pop_size,0.5*scalefactor[0],scalefactor[0],0.5*scalefactor[1],scalefactor[1])
      position_nD = [[positions[0][i], positions[1][i]] for i in range(len(positions[0]))]
    
  if mirror_neurons != None:
    # 3 lines below needed for multiple nodes # fixed !!! #
    mirror_neurons.sort(key=lambda x: x[0]) #sort by Gids and arrange related positions
    mirror_gids = [gids[0] for gids in mirror_neurons]
    mirror_pos = [pos[1] for pos in mirror_neurons]

    # add all positions together
    print('mirror neurons!   original position_nD: ',len(position_nD),'   ',position_nD[:3])
    print('mirror neurons!  mirror positions: ',len(mirror_pos),'  ',mirror_pos[:3])
    position_nD = position_nD + mirror_pos
    print('positions len for fake including mirrors: ',len(position_nD))

  if fake == 0:
    # fake == 0 is the normal case, where actual iaf neurons are instantiated
    if nucleus=='GPi_fake':
      element = 'parrot_neuron'
    else:
      nest.SetDefaults('iaf_psc_alpha_multisynapse', bg_params['common_iaf'])
      nest.SetDefaults('iaf_psc_alpha_multisynapse', bg_params[nucleus+'_iaf'])
      nest.SetDefaults('iaf_psc_alpha_multisynapse', {"I_e": bg_params['Ie'+nucleus]})
      #element = 'iaf_psc_alpha_multisynapse'
      if nucleus == 'MSN':
        nest.CopyModel('iaf_psc_alpha_multisynapse', 'msn_d1')
        nest.CopyModel('iaf_psc_alpha_multisynapse', 'msn_d2')
        element = ['msn_d1','msn_d2']
        position_nD = position_nD[:int(len(position_nD)/2)]
      else:
        element = 'iaf_psc_alpha_multisynapse'
  else:
    # when fake > 0, parrot neurons instantiated
    element = 'parrot_neuron'
  layer_gid = ntop.CreateLayer({'positions': position_nD, 'elements': element, 'extent':my_extent, 'center':my_center, 'edge_wrap': True})
  print (len(position_nD))
  if nucleus =='MSN':
    save_layers_position(nucleus, layer_gid, np.array(position_nD*2))
  else:
    save_layers_position(nucleus, layer_gid, np.array(position_nD))
  if fake > 0:
    # when fake > 0, parrot neurons are connected to poisson generators firing at `fake`Hz
    ## adding 2 lines below just in case (fix for multiple nodes)
    my_post = list(nest.GetNodes(layer_gid)[0])
    my_post.sort()

    if True:
        poisson = nest.Create('poisson_generator', 1) 
        nest.SetStatus(poisson, {'rate': fake})
        poisson_string = poisson * pop_size
        nest.Connect(pre=poisson_string, post=my_post[0:pop_size], conn_spec={'rule':'one_to_one'})
    ################################
    if mirror_neurons != None:
      print('special handling of '+ nucleus + ' input layer => the remaining neurons will be connected to the original ctx neurons')
      print('connecting mirror neurons of len: ',len(mirror_gids),' to ',nucleus)
      nest.Connect(pre=mirror_gids, post=my_post[-len(mirror_gids):], conn_spec={'rule':'one_to_one'},syn_spec={'delay':10.}) ## added delay !!!!
 
  return layer_gid

#-------------------------------------------------------------------------------
# Establishes a topological connection between two populations
# bg_params: basal ganglia parameters
# nType : a string 'ex' or 'in', defining whether it is excitatory or inhibitory
# bg_layers : dictionary of basal ganglia layers
# nameTgt, nameSrc : strings naming the populations, as defined in NUCLEI list
# projType : type of projections. For the moment: 'focused' (only channel-to-channel connection) and
#            'diffuse' (all-to-one with uniform distribution)
# redundancy, RedundancyType : contrains the inDegree - see function `connect` for details
# LCGDelays : shall we use the delays obtained by (Liénard, Cos, Girard, in prep) or not (default = True)
# gain : allows to amplify the weight normally deduced from LG14
# stochastic_delays: to enable stochasticity in the axonal delays
# spreads: a 2-item list specifying the radius of focused and diffuse projections
#-------------------------------------------------------------------------------
def connect_layers_bg(bg_params, nType, bg_layers, nameSrc, nameTgt, projType, redundancy, RedundancyType, LCGDelays=True, gain=1., stochastic_delays=None, spreads=None, verbose=False,scalefactor=[1,1]):
  def printv(text):
    if verbose:
      print(text)

  printv("\n* connecting "+nameSrc+" -> "+nameTgt+" with "+nType+" "+projType+" connection")

  recType = {'AMPA':1,'NMDA':2,'GABA':3}

  if RedundancyType == 'inDegreeAbs':
    # inDegree is already provided in the right form
    inDegree = float(redundancy)
  elif RedundancyType == 'outDegreeAbs':
    #### fractional outDegree is expressed as a fraction of max axo-dendritic contacts
    inDegree = get_frac_bg(bg_params, 1./redundancy, nameSrc, nameTgt, bg_params['count'+nameSrc], bg_params['count'+nameTgt], verbose=verbose)
  elif RedundancyType == 'outDegreeCons':
    #### fractional outDegree is expressed as a ratio of min/max axo-dendritic contacts
    inDegree = get_frac_bg(bg_params, redundancy, nameSrc, nameTgt, bg_params['count'+nameSrc], bg_params['count'+nameTgt], useMin=True, verbose=verbose)
  else:
    raise KeyError('`RedundancyType` should be one of `inDegreeAbs`, `outDegreeAbs`, or `outDegreeCons`.')

  # check if in degree acceptable (not larger than number of neurons in the source nucleus)
  if projType == 'focused' and inDegree > bg_params['nb'+nameSrc]:
    printv("/!\ WARNING: required 'in degree' ("+str(inDegree)+") larger than number of neurons in individual source channels ("+str(bg_params['nb'+nameSrc])+"), thus reduced to the latter value")
    inDegree = bg_params['nb'+nameSrc]
  if projType == 'diffuse' and inDegree > bg_params['nb'+nameSrc]:
    printv("/!\ WARNING: required 'in degree' ("+str(inDegree)+") larger than number of neurons in the overall source population ("+str(bg_params['nb'+nameSrc])+"), thus reduced to the latter value")
    inDegree = bg_params['nb'+nameSrc]

  if inDegree == 0.:
    printv("/!\ WARNING: non-existent connection strength, will skip")
    return

  global AMPASynapseCounter_bg

  # prepare receptor type lists:
  if nType == 'ex':
    lRecType = ['AMPA','NMDA']
    AMPASynapseCounter_bg = AMPASynapseCounter_bg + 1
    lbl = AMPASynapseCounter_bg # needs to add NMDA later
  elif nType == 'AMPA':
    lRecType = ['AMPA']
    lbl = 0
  elif nType == 'NMDA':
    lRecType = ['NMDA']
    lbl = 0
  elif nType == 'in':
    lRecType = ['GABA']
    lbl = 0
  else:
    raise KeyError('Undefined connexion type: '+nType)

  # compute the global weight of the connection, for each receptor type:
  if nameSrc=='GPe' and nameTgt=='STN':
    W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, bg_params['GGPe_STN'], verbose=verbose)
  else:
    if nameSrc=='STN' and nameTgt=='GPi':
      W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, bg_params['GSTN_GPi'], verbose=verbose)
    else:
      if nameSrc=='MSN' and nameTgt=='MSN':
        W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, bg_params['GMSN_MSN'], verbose=verbose)
      else:
        if (nameSrc=='MSN' and nameTgt=='GPi') or (nameSrc=='MSN' and nameTgt=='GPe'):
          if bg_params['plastic_syn']:
            W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, bg_params['GMSN_GPx'], verbose=verbose)
          else:
            W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, gain, verbose=verbose)
        else:
          W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, gain, verbose=verbose)

  
  printv("  W="+str(W)+" and inDegree="+str(inDegree))

  # determine which transmission delay to use:
  if LCGDelays:
    delay = bg_params['tau'][nameSrc+'->'+nameTgt]
  else:
    delay = 1.
  if projType == 'focused': # if projections focused, input come only from the same channel as tgtChannel
    mass_connect_bg(bg_params, bg_layers, nameSrc, nameTgt, lbl, inDegree, recType[lRecType[0]], W[lRecType[0]], delay, spread=bg_params['spread_focused'], stochastic_delays = stochastic_delays, verbose=verbose)
  elif projType == 'diffuse': # if projections diffused, input connections are shared among each possible input channel equally
    mass_connect_bg(bg_params, bg_layers, nameSrc, nameTgt, lbl, inDegree, recType[lRecType[0]], W[lRecType[0]], delay, spread=bg_params['spread_diffuse']*max(scalefactor), stochastic_delays = stochastic_delays, verbose=verbose)

  if nType == 'ex':
    # mirror the AMPA connection with similarly connected NMDA connections
    src_idx = 0
    mass_mirror_bg(bg_params,nameSrc, nameTgt,nest.GetNodes(bg_layers[nameSrc])[src_idx], lbl, recType['NMDA'], W['NMDA'], delay, stochastic_delays = stochastic_delays)

  return W


#------------------------------------------------------------------------------
# Routine to perform the fast connection using nest built-in `connect` function
# - `bg_params` is basal ganglia parameters
# - `bg_layers` is the dictionary of basal ganglia layers
# - `sourceName` & `destName` are names of two different layers
# - `synapse_label` is used to tag connections and be able to find them quickly
#   with function `mass_mirror`, that adds NMDA on top of AMPA connections
# - `inDegree`, `receptor_type`, `weight`, `delay` are Nest connection params
# - `spread` is a parameter that affects the diffusion level of the connection
#------------------------------------------------------------------------------
def mass_connect_bg(bg_params, bg_layers, sourceName, destName, synapse_label, inDegree, receptor_type, weight, delay, spread, stochastic_delays=None, verbose=False):
  def printv(text):
    if verbose:
      print(text)

  # potential initialization of stochastic delays
  if stochastic_delays != None and delay > 0:
    printv('Using stochastic delays in mass-connect')
    low = delay * 0.5
    high = delay * 1.5
    sigma = delay * stochastic_delays
    delay =  {'distribution': 'normal_clipped', 'low': low, 'high': high, 'mu': delay, 'sigma': sigma}

  ## set default synapse model with the chosen label
  nest.SetDefaults('static_synapse_lbl', {'synapse_label': synapse_label, 'receptor_type': receptor_type})

  # creation of the topological connection dict
  conndict = {'connection_type': 'convergent',
              'mask': {'circular': {'radius': spread}},
              'synapse_model': 'static_synapse_lbl', 'weights': weight, 'delays':delay,
              'allow_oversized_mask': True, 'allow_multapses': True}

  if bg_params['plastic_syn']:
     
    nest.SetDefaults('syn_d1', {'synapse_label': synapse_label, 'receptor_type': receptor_type})
    conndict_dop_d1 = {'connection_type': 'convergent',
                       'mask': {'circular': {'radius': spread}},
                       'synapse_model': 'syn_d1', 'weights': weight*bg_params['plast_gain'], 'delays':delay, #0.62 #0.65 #0.5  #0.7 #scaling down the weights for reinforcement by plasticity.
                       'allow_oversized_mask': True, 'allow_multapses': True, 'targets':{'model':'msn_d1'}}  #changed weight by weight*0.5 for dop.

    nest.SetDefaults('syn_d2', {'synapse_label': synapse_label+1000, 'receptor_type': receptor_type})  #label+1000
    conndict_dop_d2 = {'connection_type': 'convergent',
                       'mask': {'circular': {'radius': spread}},
                       'synapse_model': 'syn_d2', 'weights': weight*bg_params['plast_gain'], 'delays':delay, #0.62 #0.65 #0.5  #0.7 #scaling down the weights for reinforcement by plasticity.
                       'allow_oversized_mask': True, 'allow_multapses': True, 'targets':{'model':'msn_d2'}}  #changed weight by weight*0.5 for dop.
    
  if (sourceName=='MSN' and destName=='GPe') or (sourceName=='MSN' and destName=='GPi'):
    conndict_d1 = {'connection_type': 'convergent',
                      'mask': {'circular': {'radius': spread}},
                      'synapse_model': 'static_synapse_lbl', 'weights': weight, 'delays':delay,
                      'allow_oversized_mask': True, 'allow_multapses': True, 'sources':{'model':'msn_d1'}}
    conndict_d2 = {'connection_type': 'convergent',
                      'mask': {'circular': {'radius': spread}},
                      'synapse_model': 'static_synapse_lbl', 'weights': weight, 'delays':delay,
                      'allow_oversized_mask': True, 'allow_multapses': True, 'sources':{'model':'msn_d2'}}

  if sourceName=='MSN' and destName=='MSN':
    conndict_MSN_d2_d2 = {'connection_type': 'convergent',
                      'mask': {'circular': {'radius': spread}},
                      'synapse_model': 'static_synapse_lbl', 'weights': weight*bg_params['syn_asymm'], 'delays':delay,
                      'allow_oversized_mask': True, 'allow_multapses': True, 'sources':{'model':'msn_d2'},'targets':{'model':'msn_d2'}}
    conndict_MSN_d2_d1 = {'connection_type': 'convergent',
                      'mask': {'circular': {'radius': spread}},
                      'synapse_model': 'static_synapse_lbl', 'weights': weight*bg_params['syn_asymm'], 'delays':delay,
                      'allow_oversized_mask': True, 'allow_multapses': True, 'sources':{'model':'msn_d2'},'targets':{'model':'msn_d1'}}
    conndict_MSN_d1_d1 = {'connection_type': 'convergent',
                      'mask': {'circular': {'radius': spread}},
                      'synapse_model': 'static_synapse_lbl', 'weights': weight, 'delays':delay,
                      'allow_oversized_mask': True, 'allow_multapses': True, 'sources':{'model':'msn_d1'},'targets':{'model':'msn_d1'}}
    conndict_MSN_d1_d2 = {'connection_type': 'convergent',
                      'mask': {'circular': {'radius': spread}},
                      'synapse_model': 'static_synapse_lbl', 'weights': weight, 'delays':delay,
                      'allow_oversized_mask': True, 'allow_multapses': True, 'sources':{'model':'msn_d1'},'targets':{'model':'msn_d2'}}
  

  # The first call ensures that all neurons in `destName`
  # have at least `int(inDegree)` incoming connections
  integer_inDegree = np.floor(inDegree)
  if integer_inDegree>0:
    printv('Adding '+str(int(integer_inDegree*bg_params['nb'+destName]))+' connections with rule `fixed_indegree`')
    integer_conndict = conndict.copy()
    integer_conndict.update({'number_of_connections': int(integer_inDegree)})
    if bg_params['plastic_syn']:
      
      integer_conndict_dop_d1 = conndict_dop_d1.copy()
      integer_conndict_dop_d1.update({'number_of_connections': int(integer_inDegree)})

      integer_conndict_dop_d2 = conndict_dop_d2.copy()
      integer_conndict_dop_d2.update({'number_of_connections': int(integer_inDegree)})
      
    
    if (sourceName=='MSN' and destName=='GPe') or (sourceName=='MSN' and destName=='GPi'):
      
      if destName=='GPi':
        integer_conndict_d1 = conndict_d1.copy()
        integer_conndict_d1.update({'number_of_connections': int(integer_inDegree*(1.-bg_params['overlap_d1d2']))})
        integer_conndict_d2 = conndict_d2.copy()
        integer_conndict_d2.update({'number_of_connections': int(integer_inDegree*bg_params['overlap_d1d2'])})
      if destName=='GPe':
        integer_conndict_d1 = conndict_d1.copy()
        integer_conndict_d1.update({'number_of_connections': int(integer_inDegree*bg_params['overlap_d1d2'])})
        integer_conndict_d2 = conndict_d2.copy()
        integer_conndict_d2.update({'number_of_connections': int(integer_inDegree*(1.-bg_params['overlap_d1d2']))})

    if (sourceName=='MSN' and destName=='MSN'):
      integer_conndict_MSN_d2_d2 = conndict_MSN_d2_d2.copy()
      integer_conndict_MSN_d2_d2.update({'number_of_connections': int(integer_inDegree*bg_params['asymmetry_1'])}) #MSN_d2 receives 50% of connctions from MSN d2 

      integer_conndict_MSN_d2_d1 = conndict_MSN_d2_d1.copy()
      integer_conndict_MSN_d2_d1.update({'number_of_connections':int(integer_inDegree*bg_params['asymmetry_1'])}) # MSN_d1 receives 50% of connctions from MSN d2 
      
      integer_conndict_MSN_d1_d1 = conndict_MSN_d1_d1.copy()
      integer_conndict_MSN_d1_d1.update({'number_of_connections':int(integer_inDegree*bg_params['asymmetry_1'])}) # MSN_d1 receives 50% of connctions from MSN d1 
      
      integer_conndict_MSN_d1_d2 = conndict_MSN_d1_d2.copy()
      integer_conndict_MSN_d1_d2.update({'number_of_connections':int(integer_inDegree*bg_params['asymmetry_2'])}) #MSN_d1 to MSN_d2 receives a small fraction of connections (asymmetry)


    if (sourceName=='CSN' and destName=='MSN') or (sourceName=='PTN' and destName=='MSN'):
      if bg_params['plastic_syn']: 
        ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], integer_conndict_dop_d1)
        ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], integer_conndict_dop_d2)

        with open('./log/weights_d1.txt','a') as file:
          file.write(sourceName+'->'+destName+':  '+str(weight * bg_params['plast_gain'])+' \n')
          file.write('at mass connect setting for stdp_dopamine_synapse_lbl d1: '+str(nest.GetDefaults('syn_d1'))+' \n')
        with open('./log/weights_d2.txt','a') as file:
          file.write(sourceName+'->'+destName+':  '+str(weight * bg_params['plast_gain'])+' \n')
          file.write('at mass connect setting for stdp_dopamine_synapse_lbl d2: '+str(nest.GetDefaults('syn_d2'))+' \n')
      else:
        ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], integer_conndict)
    else:
      if (sourceName=='MSN' and destName=='GPe') or (sourceName=='MSN' and destName=='GPi'):
        ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], integer_conndict_d1)
        ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], integer_conndict_d2)
      else:
        if (sourceName=='MSN' and destName=='MSN'):
          ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], integer_conndict_MSN_d2_d2)
          ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], integer_conndict_MSN_d2_d1)
          ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], integer_conndict_MSN_d1_d1)
          ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], integer_conndict_MSN_d1_d2)
        else:
          ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], integer_conndict)

  # The second call distributes the approximate number of remaining axonal
  # contacts at random (i.e. the remaining fractional part after the first step)
  # Why "approximate"? Because with pynest layers, there are only two ways to specify
  # the number of axons in a connection:
  #    1) with an integer, specified with respect to each source (alt. target) neurons
  #    2) as a probability
  # Here, we have a fractional part - not an integer number - so that leaves us option 2.
  # However, because the new axonal contacts are drawn at random, we will not have the
  # exact number of connections
  float_inDegree = inDegree - integer_inDegree
  remaining_connections = np.round(float_inDegree * bg_params['nb'+destName])
  if remaining_connections > 0:
    printv('Adding '+str(remaining_connections)+' remaining connections with rule `fixed_total_number`')
    float_conndict = conndict.copy()
    float_conndict.update({'kernel': 1. / (bg_params['nb'+sourceName] * float(remaining_connections))})
    ntop.ConnectLayers(bg_layers[sourceName], bg_layers[destName], float_conndict)

#------------------------------------------------------------------------------
# Routine to duplicate a connection made with a specific receptor, with another
# receptor (typically to add NMDA connections to existing AMPA connections)
# - `source` & `synapse_label` should uniquely define the connections of
#   interest - typically, they are the same as in the call to `mass_connect`
# - `receptor_type`, `weight`, `delay` are Nest connection params
#------------------------------------------------------------------------------
def mass_mirror_bg(bg_params, nameSrc, nameTgt,source, synapse_label, receptor_type, weight, delay, stochastic_delays, verbose=False):
  def printv(text):
    if verbose:
      print(text)

  # find all AMPA connections for the given projection type
  printv('looking for AMPA connections to mirror with NMDA...\n')
  ampa_conns = nest.GetConnections(source=source, synapse_label=synapse_label)
  # in rare cases, there may be no connections, guard against that
  if ampa_conns:
    # extract just source and target GID lists, all other information is irrelevant here
    printv('found '+str(len(ampa_conns))+' AMPA connections\n')
    if stochastic_delays != None and delay > 0:
      printv('Using stochastic delays in mass-miror')
      delay = np.array(nest.GetStatus(ampa_conns, keys=['delay'])).flatten()
    src, tgt, _, _, _ = zip(*ampa_conns)
    #### adding plastic synapses ########
    #####################################
    if (nameSrc=='CSN' and nameTgt=='MSN') or (nameSrc=='PTN' and nameTgt=='MSN'):
      
      if bg_params['plastic_syn']:
        nest.Connect(src, tgt, 'one_to_one',{'model': 'syn_d1',
                    'receptor_type': receptor_type, 'weight': weight*bg_params['plast_gain'], 'delay':delay})
        
        ampa_conns_d2 = nest.GetConnections(source=source, synapse_label=synapse_label+1000)
        if ampa_conns_d2:
          src_d2, tgt_d2, _, _, _ = zip(*ampa_conns_d2)
          nest.Connect(src_d2, tgt_d2, 'one_to_one',{'model': 'syn_d2',
                      'receptor_type': receptor_type, 'weight': weight*bg_params['plast_gain'], 'delay':delay})

        with open('./log/weights_d1.txt','a') as file:
          file.write('mirror NMDA d1, weight: '+str(weight*bg_params['plast_gain'])+' \n')
          file.write('at mirror connect stdp_dopamine_synapse_lbl d1: '+str(nest.GetDefaults('syn_d1'))+' \n')
        with open('./log/weights_d2.txt','a') as file:
          file.write('mirror NMDA d2, weight: '+str(weight*bg_params['plast_gain'])+' \n')
          file.write('at mirror connect stdp_dopamine_synapse_lbl d2: '+str(nest.GetDefaults('syn_d2'))+' \n')
        
      else:
        nest.Connect(src, tgt, 'one_to_one',
                 {'model': 'static_synapse_lbl',
                  'synapse_label': synapse_label, # tag with the same number (doesn't matter)
                  'receptor_type': receptor_type, 'weight': weight, 'delay':delay})
    #########################################
    else:
      nest.Connect(src, tgt, 'one_to_one',
                 {'model': 'static_synapse_lbl',
                  'synapse_label': synapse_label, # tag with the same number (doesn't matter)
                  'receptor_type': receptor_type, 'weight': weight, 'delay':delay})

  
#-------------------------------------------------------------------------------
# Helper function to set a basal ganglia internal projection
# computes the inDegree as a fraction of maximal possible inDegree
# `FractionalOutDegree` is the outDegree, expressed as a fraction
#-------------------------------------------------------------------------------
def get_frac_bg(bg_params, FractionalOutDegree, nameSrc, nameTgt, cntSrc, cntTgt, useMin=False, verbose=False):
  if useMin == False:
    # 'FractionalOutDegree' is taken to be relative to the maximal number of axo-dendritic contacts
    inDegree = get_input_range_bg(bg_params, nameSrc, nameTgt, cntSrc, cntTgt, verbose=verbose)[1] * FractionalOutDegree
  else:
    # 'FractionalOutDegree' is taken to be relative to the maximal number of axo-dendritic contacts and their minimal number
    r = get_input_range_bg(bg_params, nameSrc, nameTgt, cntSrc, cntTgt, verbose=verbose)
    inDegree = (r[1] - r[0]) * FractionalOutDegree + r[0]
  if verbose:
    print('\tConverting the fractional outDegree of '+nameSrc+' -> '+nameTgt+' from '+str(FractionalOutDegree)+' to inDegree neuron count: '+str(round(inDegree, 2))+' (relative to minimal value possible? '+str(useMin)+')')
  return inDegree

#-------------------------------------------------------------------------------
# Helper function to set a basal ganglia internal projection
# computes the weight of a connection, based on LG14 parameters
#-------------------------------------------------------------------------------
def computeW_bg(bg_params, listRecType, nameSrc, nameTgt, inDegree, gain=1.,verbose=False):
  recType = {'AMPA':1,'NMDA':2,'GABA':3}
  nu = get_input_range_bg(bg_params, nameSrc, nameTgt, bg_params['count'+nameSrc], bg_params['count'+nameTgt], verbose=verbose)[1]
  if verbose:
    print('\tCompare with the effective chosen inDegree   : '+str(inDegree))

  # attenuation due to the distance from the receptors to the soma of tgt:
  LX=bg_params['lx'][nameTgt]*np.sqrt((4.*bg_params['Ri'])/(bg_params['dx'][nameTgt]*bg_params['Rm']))
  attenuation = np.cosh(LX*(1-bg_params['distcontact'][nameSrc+'->'+nameTgt])) / np.cosh(LX)

  w={}
  for r in listRecType:
    w[r] = nu / float(inDegree) * attenuation * bg_params['wPSP'][recType[r]-1] * gain
  return w

#-------------------------------------------------------------------------------
# Helper function to set a basal ganglia internal projection
# returns the minimal & maximal numbers of distinct input neurons for one connection
#-------------------------------------------------------------------------------
def get_input_range_bg(bg_params, nameSrc, nameTgt, cntSrc, cntTgt, verbose=False):
  if nameSrc=='CSN' or nameSrc=='PTN':
    nu = bg_params['alpha'][nameSrc+'->'+nameTgt]
    nu0 = 0
    if verbose:
      print('\tMaximal number of distinct input neurons (nu): '+str(nu))
      print('\tMinimal number of distinct input neurons     : unknown (set to 0)')
  else:
    nu = cntSrc / float(cntTgt) * bg_params['ProjPercent'][nameSrc+'->'+nameTgt] * bg_params['alpha'][nameSrc+'->'+nameTgt]
    nu0 = cntSrc / float(cntTgt) * bg_params['ProjPercent'][nameSrc+'->'+nameTgt]
    if verbose:
      print('\tMaximal number of distinct input neurons (nu): '+str(nu))
      print('\tMinimal number of distinct input neurons     : '+str(nu0))
  return [nu0, nu]


def connect_GPi2d_GPi3d(GPi2d,GPi3d): #connect GPi layer to fake GPi
  my_pre = nest.GetNodes(GPi2d)[0]
  my_post = nest.GetNodes(GPi3d)[0]
  nest.SetDefaults('static_synapse',{'receptor_type':0})
  nest.Connect(pre=my_pre, post=my_post, conn_spec={'rule':'one_to_one'})
