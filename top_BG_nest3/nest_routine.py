#!/usr/bin/env python
# -*- coding: utf-8 -*-

##
## nest_routine.py  —  NEST 3 compatible
##
## Changes from NEST 2:
##   - nest.topology (ntop) removed; spatial layers use nest.spatial.free + nest.Create
##   - spike_detector  →  spike_recorder  (params: record_to='ascii')
##   - nest.GetNodes(layer)[0]  →  layer  (NodeCollection directly)
##   - nest.GetKernelStatus(['key'])[0]  →  nest.GetKernelStatus('key')
##   - grng_seed + rng_seeds  →  single rng_seed
##   - ntop.ConnectLayers(src, tgt, conndict)  →  nest.Connect(src, tgt, conn_spec, syn_spec)
##   - connection_type/number_of_connections  →  rule='fixed_indegree', indegree=N
##   - kernel (probability)  →  rule='pairwise_bernoulli', p=...
##   - sources/targets model filters  →  NodeCollection slicing (MSN[:half] / MSN[half:])
##   - SynapseCollection unpacking via zip(*conns)  →  conns.get([...])
##   - stochastic delays dict  →  nest.math.redraw(nest.random.normal(...))
##   - SetStatus([psg], ...)  →  psg.set(...)
##

import nest
import numpy as np
import math
import time
import collections
import os
import random
import itertools

pyrngs = []


# -------------------------------------------------------------------------------
# NEST initialization
# -------------------------------------------------------------------------------
def initialize_nest(sim_params):
    nest.set_verbosity("M_WARNING")
    nest.SetKernelStatus({"overwrite_files": sim_params['overwrite_files']})
    nest.SetKernelStatus({'local_num_threads': int(sim_params['nbcpu'])})
    nest.SetKernelStatus({"data_path": 'log'})
    if sim_params['dt'] != '0.1':
        nest.SetKernelStatus({'resolution': float(sim_params['dt'])})
    # NEST 3: single rng_seed replaces grng_seed + rng_seeds
    N_vp = nest.GetKernelStatus('total_num_virtual_procs')
    print('Number of virtual processes: ', N_vp)
    global pyrngs
    pyrngs = [np.random.RandomState(s) for s in range(sim_params['msd'], sim_params['msd'] + N_vp)]
    nest.SetKernelStatus({'rng_seed': sim_params['msd']})


# -------------------------------------------------------------------------------
# Spike recorder (NEST 3: spike_detector → spike_recorder)
# -------------------------------------------------------------------------------
def layer_spike_detector(layer_gid, layer_name, ignore_time, params=None):
    if params is None:
        params = {}
    p = dict(params)
    p.update({'record_to': 'ascii', 'label': layer_name, 'start': float(ignore_time)})
    print('spike recorder for ' + layer_name)
    detector = nest.Create("spike_recorder", params=p)
    nest.Connect(layer_gid, detector)
    return detector


def get_events_from_ascii(detector):
    filenames = detector.get('filenames')
    if isinstance(filenames, str):
        filenames = [filenames]
    senders, times = [], []
    for fname in filenames:
        if fname and os.path.isfile(fname):
            with open(fname) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            senders.append(int(parts[0]))
                            times.append(float(parts[1]))
                        except ValueError:
                            pass
    return {'senders': np.array(senders, dtype=int), 'times': np.array(times)}


def average_fr(detector, simDuration, n):
    # NEST 3: use NodeCollection .get() API
    return detector.get('n_events') / (float(simDuration) * float(n) / 1000.)


def instantaneous_fr(detector, sim_start, sim_end, n, dt=0.1):
    At = []
    filenames = detector.get('filenames')
    if isinstance(filenames, str):
        filenames = [filenames]
    spike_times = []
    for fname in filenames:
        if fname and os.path.isfile(fname):
            with open(fname) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            spike_times.append(float(parts[1]))
                        except ValueError:
                            pass  # skip header row
    spike_times = np.array(spike_times)
    for i in np.arange(int(sim_start), int(sim_end), dt):
        a = ((i <= spike_times) & (spike_times < (i + 1))).sum()
        At.append(a)
    At = np.array(At) / float(n) * 1000.
    return At


# -------------------------------------------------------------------------------
# NEST 3: layer_gid is a NodeCollection — len() replaces GetNodes()[0]
# -------------------------------------------------------------------------------
def count_layer(layer_gid):
    return len(layer_gid)


def get_connection(gids):
    return nest.GetConnections(gids)


def get_columns_data(layer_name, circle_center, radius_small, my_area=None):
    gid_pos = np.loadtxt('./log/' + layer_name + '.txt')
    print('gids and pos  ', len(gid_pos))
    circle_gids = []
    for i in np.arange(len(circle_center)):
        idx = np.where(np.linalg.norm([(gid_pos[:, 1] - circle_center[i][0]),
                                       (gid_pos[:, 2] - circle_center[i][1])], axis=0) <= radius_small)[0]
        print('number of neurons in channel ', str(i), 'for ', layer_name, ' : ', str(len(idx)))
        circle_gids.append([[int(x[0]), x[1:].tolist()] for x in gid_pos[idx, :]])
    return circle_gids


AMPASynapseCounter_bg = 0


def hex_corner(center, size, i):
    angle_deg = 60 * i - 30
    angle_rad = np.pi / 180 * angle_deg
    return [center[0] + size * np.cos(angle_rad), center[1] + size * np.sin(angle_rad)]


def get_channel_centers(bg_params, hex_center=[0, 0], ci=6, hex_radius=0.240):
    center_aux = []
    if bg_params['channels']:
        if len(bg_params['circle_center']) == 0:
            for i in np.arange(ci):
                x_y = hex_corner(hex_center, hex_radius, i)
                center_aux.append(x_y)
            np.savetxt('./log/centers.txt', center_aux)
            print('generated centers: ', center_aux)
    return center_aux


def apply_fixed_stimulus(time_start, time_stop, psg, stimulus):
    # NEST 3: NodeCollection.set() replaces nest.SetStatus([node], ...)
    psg.set({'rate': stimulus, 'start': time_start, 'stop': time_stop})


def get_targets_mean_rates(columns_gids, spike_detector, idx_dt, t_start, t_end):
    targets_rates = []
    for col in np.arange(5):
        gid_l = [gid[0] for gid in columns_gids[col]]
        mask = np.isin(spike_detector['senders'][idx_dt], gid_l)
        mean_act = (len(spike_detector['times'][idx_dt][mask])
                    / ((t_end - t_start) * float(len(columns_gids[col]))) * 1000.)
        print('mean firing rate for col ', col, ':  ', mean_act)
        targets_rates.append(mean_act)
    return targets_rates


def create_psg_channels(syn_weight, syn_delay, channels_nb):
    psg = nest.Create('poisson_generator', channels_nb)
    syn = {'weight': syn_weight, 'delay': syn_delay}
    return psg, syn


def connect_psg_to_channels(columns_gids, psg, syn):
    circle_j_gids_nb = {}
    for j in np.arange(len(columns_gids)):
        circle_j_gids = sorted([k[0] for k in columns_gids[j]])
        print('circle  ', j, ' contains #: ', len(circle_j_gids))
        circle_j_gids_nb[str(j)] = len(circle_j_gids)
        # NEST 3: psg[j] is NodeCollection; target GIDs wrapped in NodeCollection
        nest.Connect(psg[j], nest.NodeCollection(circle_j_gids), syn_spec=syn)


def grid_positions(nbCh, sim_pts, a0, a1, b0, b1):
    n = int(sim_pts * nbCh)
    n_squared = np.ceil(np.sqrt(n))
    coord = [[x / n_squared * a1 - a0, y / n_squared * b1 - b0]
              for x in np.arange(0, n_squared, dtype=float)
              for y in np.arange(0, n_squared, dtype=float)]
    if len(coord) > n:
        coord = np.array(coord)[np.sort(pyrngs[0].choice(range(len(coord)), size=n, replace=False))].tolist()
    aux_x = [coord[i][0] for i in range(len(coord))]
    aux_y = [coord[i][1] for i in range(len(coord))]
    return [aux_x, aux_y]


def grid_uniform_positions(number_of_neurons):
    pos = [[np.random.uniform(-0.5, 0.5), np.random.uniform(-0.5, 0.5)] for j in range(number_of_neurons)]
    return pos


def save_layers_position(layer_name, layer_gid, positions):
    # NEST 3: NodeCollection.tolist() returns list of node IDs
    node_ids = np.array(layer_gid.tolist())
    gid_and_positions = np.column_stack((node_ids, positions))
    np.savetxt('log/' + layer_name + '.txt', gid_and_positions, fmt='%1.3f')
    if layer_name == 'MSN':
        half = int(len(gid_and_positions) / 2)
        np.savetxt('log/' + layer_name + '_d1.txt', gid_and_positions[:half], fmt='%1.3f')
        np.savetxt('log/' + layer_name + '_d2.txt', gid_and_positions[half:], fmt='%1.3f')


# -------------------------------------------------------------------------------
# Create a topological layer
# NEST 3: ntop.CreateLayer → nest.Create with nest.spatial.free positions
# MSN special case: d1 and d2 created as separate NodeCollections then merged
# -------------------------------------------------------------------------------
def create_layers_bg(bg_params, nucleus, fake=0, mirror_neurons=None, mirror_pos=None, scalefactor=[1, 1]):

    my_extent = [1. * int(scalefactor[0]) + 1., 1. * int(scalefactor[1]) + 1.]

    if mirror_neurons is None:
        if nucleus == 'GPi_fake':
            pop_size = int(bg_params['nb' + nucleus[:3]])
        else:
            pop_size = int(bg_params['nb' + nucleus])
    else:
        pop_size = int(bg_params['nb' + nucleus]) - len(mirror_neurons)

    print('population size for ' + nucleus + ': ' + str(pop_size))

    if nucleus == 'GPi_fake':
        positions_z = pyrngs[0].uniform(0., 0.5, pop_size).tolist()
        positions = np.loadtxt('./log/' + nucleus[:3] + '.txt')
        position_nD = [[positions[i][1], positions[i][2], positions_z[i]] for i in range(len(positions))]
        my_extent = my_extent + [1.]
        print('positions GPi_fake: ', position_nD[:10])
    else:
        if nucleus == 'MSN':
            positions = grid_uniform_positions(pop_size)
            position_nD = positions
        else:
            if nucleus == 'GPi' or nucleus == 'STN':
                positions = grid_positions(1, pop_size,
                                           0.4 * scalefactor[0], scalefactor[0] - 0.1,
                                           0.4 * scalefactor[1], scalefactor[1] - 0.1)
            else:
                positions = grid_positions(1, pop_size,
                                           0.5 * scalefactor[0], scalefactor[0],
                                           0.5 * scalefactor[1], scalefactor[1])
            position_nD = [[positions[0][i], positions[1][i]] for i in range(len(positions[0]))]

    if mirror_neurons is not None:
        mirror_neurons.sort(key=lambda x: x[0])
        mirror_gids = [gids[0] for gids in mirror_neurons]
        mirror_pos = [pos[1] for pos in mirror_neurons]
        print('mirror neurons!   original position_nD: ', len(position_nD), '   ', position_nD[:3])
        print('mirror neurons!  mirror positions: ', len(mirror_pos), '  ', mirror_pos[:3])
        position_nD = position_nD + mirror_pos
        print('positions len for fake including mirrors: ', len(position_nD))

    if fake == 0:
        if nucleus == 'GPi_fake':
            spatial_pos = nest.spatial.free(position_nD, extent=my_extent, edge_wrap=True)
            layer_gid = nest.Create('parrot_neuron', positions=spatial_pos)

        elif nucleus == 'MSN':
            # NEST 3: d1 and d2 kept separate — NodeCollections with different models
            # cannot be joined with + in NEST 3 even when positions match.
            nest.SetDefaults('iaf_psc_alpha_multisynapse', bg_params['common_iaf'])
            nest.SetDefaults('iaf_psc_alpha_multisynapse', bg_params[nucleus + '_iaf'])
            nest.SetDefaults('iaf_psc_alpha_multisynapse', {"I_e": bg_params['Ie' + nucleus]})
            nest.CopyModel('iaf_psc_alpha_multisynapse', 'msn_d1')
            nest.CopyModel('iaf_psc_alpha_multisynapse', 'msn_d2')
            half = int(pop_size / 2)
            pos_half = position_nD[:half]
            spatial_d1 = nest.spatial.free(pos_half, extent=my_extent, edge_wrap=True)
            spatial_d2 = nest.spatial.free(pos_half, extent=my_extent, edge_wrap=True)
            layer_d1 = nest.Create('msn_d1', positions=spatial_d1)
            layer_d2 = nest.Create('msn_d2', positions=spatial_d2)
            save_layers_position('MSN_d1', layer_d1, np.array(pos_half))
            save_layers_position('MSN_d2', layer_d2, np.array(pos_half))
            d1_data = np.column_stack((np.array(layer_d1.tolist()), np.array(pos_half)))
            d2_data = np.column_stack((np.array(layer_d2.tolist()), np.array(pos_half)))
            np.savetxt('log/MSN.txt', np.vstack((d1_data, d2_data)), fmt='%1.3f')
            return layer_d1, layer_d2

        else:
            nest.SetDefaults('iaf_psc_alpha_multisynapse', bg_params['common_iaf'])
            nest.SetDefaults('iaf_psc_alpha_multisynapse', bg_params[nucleus + '_iaf'])
            nest.SetDefaults('iaf_psc_alpha_multisynapse', {"I_e": bg_params['Ie' + nucleus]})
            spatial_pos = nest.spatial.free(position_nD, extent=my_extent, edge_wrap=True)
            layer_gid = nest.Create('iaf_psc_alpha_multisynapse', positions=spatial_pos)
    else:
        spatial_pos = nest.spatial.free(position_nD, extent=my_extent, edge_wrap=True)
        layer_gid = nest.Create('parrot_neuron', positions=spatial_pos)

    print(len(position_nD))
    save_layers_position(nucleus, layer_gid, np.array(position_nD))

    if fake > 0:
        # NEST 3: one poisson_generator broadcast to all parrot neurons (1 source → all_to_all = broadcast)
        poisson = nest.Create('poisson_generator', 1)
        poisson.set({'rate': fake})
        nest.Connect(poisson, layer_gid, conn_spec={'rule': 'all_to_all'})

        if mirror_neurons is not None:
            print('special handling of ' + nucleus + ' => connecting mirror neurons')
            print('connecting mirror neurons of len: ', len(mirror_gids), ' to ', nucleus)
            nest.Connect(nest.NodeCollection(mirror_gids), layer_gid[-len(mirror_gids):],
                         conn_spec={'rule': 'one_to_one'}, syn_spec={'delay': 10.})

    return layer_gid


# -------------------------------------------------------------------------------
# Topological connection between two populations
# NEST 3: conn_spec + syn_spec replace the NEST 2 conndict passed to ConnectLayers
# -------------------------------------------------------------------------------
def connect_layers_bg(bg_params, nType, bg_layers, nameSrc, nameTgt, projType, redundancy,
                      RedundancyType, LCGDelays=True, gain=1., stochastic_delays=None,
                      spreads=None, verbose=False, scalefactor=[1, 1]):
    def printv(text):
        if verbose:
            print(text)

    printv("\n* connecting " + nameSrc + " -> " + nameTgt + " with " + nType + " " + projType + " connection")

    recType = {'AMPA': 1, 'NMDA': 2, 'GABA': 3}

    if RedundancyType == 'inDegreeAbs':
        inDegree = float(redundancy)
    elif RedundancyType == 'outDegreeAbs':
        inDegree = get_frac_bg(bg_params, 1. / redundancy, nameSrc, nameTgt,
                                bg_params['count' + nameSrc], bg_params['count' + nameTgt], verbose=verbose)
    elif RedundancyType == 'outDegreeCons':
        inDegree = get_frac_bg(bg_params, redundancy, nameSrc, nameTgt,
                                bg_params['count' + nameSrc], bg_params['count' + nameTgt],
                                useMin=True, verbose=verbose)
    else:
        raise KeyError('`RedundancyType` should be one of `inDegreeAbs`, `outDegreeAbs`, or `outDegreeCons`.')

    if projType == 'focused' and inDegree > bg_params['nb' + nameSrc]:
        printv("/!\\ WARNING: required 'in degree' larger than source population, reduced")
        inDegree = bg_params['nb' + nameSrc]
    if projType == 'diffuse' and inDegree > bg_params['nb' + nameSrc]:
        printv("/!\\ WARNING: required 'in degree' larger than source population, reduced")
        inDegree = bg_params['nb' + nameSrc]

    if inDegree == 0.:
        printv("/!\\ WARNING: non-existent connection strength, will skip")
        return

    global AMPASynapseCounter_bg

    if nType == 'ex':
        lRecType = ['AMPA', 'NMDA']
        AMPASynapseCounter_bg += 1
        lbl = AMPASynapseCounter_bg
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
        raise KeyError('Undefined connexion type: ' + nType)

    if nameSrc == 'GPe' and nameTgt == 'STN':
        W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, bg_params['GGPe_STN'], verbose=verbose)
    elif nameSrc == 'STN' and nameTgt == 'GPi':
        W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, bg_params['GSTN_GPi'], verbose=verbose)
    elif nameSrc == 'MSN' and nameTgt == 'MSN':
        W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, bg_params['GMSN_MSN'], verbose=verbose)
    elif (nameSrc == 'MSN' and nameTgt == 'GPi') or (nameSrc == 'MSN' and nameTgt == 'GPe'):
        if bg_params['plastic_syn']:
            W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, bg_params['GMSN_GPx'], verbose=verbose)
        else:
            W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, gain, verbose=verbose)
    else:
        W = computeW_bg(bg_params, lRecType, nameSrc, nameTgt, inDegree, gain, verbose=verbose)

    printv("  W=" + str(W) + " and inDegree=" + str(inDegree))

    if LCGDelays:
        delay = bg_params['tau'][nameSrc + '->' + nameTgt]
    else:
        delay = 1.

    spread = (bg_params['spread_focused'] if projType == 'focused'
              else bg_params['spread_diffuse'] * max(scalefactor))

    # NEST 3: use CollocatedSynapses to co-localise AMPA+NMDA on the same pre→post pairs
    mass_connect_bg(bg_params, bg_layers, nameSrc, nameTgt, lbl, inDegree,
                    recType[lRecType[0]], W[lRecType[0]], delay,
                    spread=spread, stochastic_delays=stochastic_delays, verbose=verbose,
                    nmda_receptor_type=(recType['NMDA'] if nType == 'ex' else None),
                    nmda_weight=(W['NMDA'] if nType == 'ex' else None))

    return W


# -------------------------------------------------------------------------------
# Fast spatial connection
# NEST 3: ntop.ConnectLayers → nest.Connect with conn_spec + syn_spec
#         source/target model filters → NodeCollection slicing (MSN[:half] / MSN[half:])
# -------------------------------------------------------------------------------
def mass_connect_bg(bg_params, bg_layers, sourceName, destName, synapse_label, inDegree,
                    receptor_type, weight, delay, spread, stochastic_delays=None, verbose=False,
                    nmda_receptor_type=None, nmda_weight=None):
    """Connect two populations.
    When nmda_receptor_type/nmda_weight are provided, AMPA and NMDA synapses are
    co-localised on the same pre→post pairs using nest.CollocatedSynapses.
    """
    def printv(text):
        if verbose:
            print(text)

    # NEST 3 stochastic delays via nest.math parameter objects
    if stochastic_delays is not None and delay > 0:
        printv('Using stochastic delays in mass-connect')
        sigma = delay * stochastic_delays
        low = delay * 0.5
        high = delay * 1.5
        delay_param = nest.math.redraw(nest.random.normal(mean=delay, std=sigma), min=low, max=high)
    else:
        delay_param = delay

    # Build syn_spec: CollocatedSynapses for excitatory (AMPA+NMDA), plain dict for inhibitory
    ampa_spec = {'synapse_model': 'static_synapse_lbl', 'synapse_label': synapse_label,
                 'receptor_type': receptor_type, 'weight': weight, 'delay': delay_param}
    if nmda_receptor_type is not None:
        nmda_spec = {'synapse_model': 'static_synapse_lbl', 'synapse_label': synapse_label,
                     'receptor_type': nmda_receptor_type, 'weight': nmda_weight, 'delay': delay_param}
        base_syn_spec = nest.CollocatedSynapses(ampa_spec, nmda_spec)
    else:
        base_syn_spec = ampa_spec

    integer_inDegree = int(np.floor(inDegree))

    if integer_inDegree > 0:
        printv('Adding ' + str(integer_inDegree * bg_params['nb' + destName]) + ' connections (fixed_indegree)')

        base_conn_spec = {
            'rule': 'fixed_indegree',
            'indegree': integer_inDegree,
            'mask': {'circular': {'radius': spread}},
            'allow_oversized_mask': True,
            'allow_multapses': True,
        }

        if (sourceName == 'CSN' and destName == 'MSN') or (sourceName == 'PTN' and destName == 'MSN'):
            if bg_params['plastic_syn']:
                w_plast = weight * bg_params['plast_gain']
                ampa_d1 = {'synapse_model': 'syn_d1', 'synapse_label': synapse_label,
                           'receptor_type': receptor_type, 'weight': w_plast, 'delay': delay_param}
                ampa_d2 = {'synapse_model': 'syn_d2', 'synapse_label': synapse_label + 1000,
                           'receptor_type': receptor_type, 'weight': w_plast, 'delay': delay_param}
                if nmda_receptor_type is not None:
                    w_plast_nmda = nmda_weight * bg_params['plast_gain']
                    nmda_d1 = {'synapse_model': 'syn_d1', 'synapse_label': synapse_label,
                               'receptor_type': nmda_receptor_type, 'weight': w_plast_nmda, 'delay': delay_param}
                    nmda_d2 = {'synapse_model': 'syn_d2', 'synapse_label': synapse_label + 1000,
                               'receptor_type': nmda_receptor_type, 'weight': w_plast_nmda, 'delay': delay_param}
                    syn_d1 = nest.CollocatedSynapses(ampa_d1, nmda_d1)
                    syn_d2 = nest.CollocatedSynapses(ampa_d2, nmda_d2)
                else:
                    syn_d1 = ampa_d1
                    syn_d2 = ampa_d2
                nest.Connect(bg_layers[sourceName], bg_layers['MSN_d1'],
                             conn_spec=base_conn_spec, syn_spec=syn_d1)
                nest.Connect(bg_layers[sourceName], bg_layers['MSN_d2'],
                             conn_spec=base_conn_spec, syn_spec=syn_d2)
                with open('./log/weights_d1.txt', 'a') as f:
                    f.write(sourceName + '->' + destName + ':  ' + str(w_plast) + ' \n')
                with open('./log/weights_d2.txt', 'a') as f:
                    f.write(sourceName + '->' + destName + ':  ' + str(w_plast) + ' \n')
            else:
                nest.Connect(bg_layers[sourceName], bg_layers['MSN_d1'],
                             conn_spec=base_conn_spec, syn_spec=base_syn_spec)
                nest.Connect(bg_layers[sourceName], bg_layers['MSN_d2'],
                             conn_spec=base_conn_spec, syn_spec=base_syn_spec)

        elif (sourceName == 'MSN' and destName == 'GPe') or (sourceName == 'MSN' and destName == 'GPi'):
            if destName == 'GPi':
                n_d1 = int(integer_inDegree * (1. - bg_params['overlap_d1d2']))
                n_d2 = int(integer_inDegree * bg_params['overlap_d1d2'])
            else:  # GPe
                n_d1 = int(integer_inDegree * bg_params['overlap_d1d2'])
                n_d2 = int(integer_inDegree * (1. - bg_params['overlap_d1d2']))
            nest.Connect(bg_layers['MSN_d1'], bg_layers[destName],
                         conn_spec=dict(base_conn_spec, indegree=n_d1), syn_spec=base_syn_spec)
            nest.Connect(bg_layers['MSN_d2'], bg_layers[destName],
                         conn_spec=dict(base_conn_spec, indegree=n_d2), syn_spec=base_syn_spec)

        elif sourceName == 'MSN' and destName == 'MSN':
            n_a1 = int(integer_inDegree * bg_params['asymmetry_1'])
            n_a2 = int(integer_inDegree * bg_params['asymmetry_2'])
            asym_syn = dict(ampa_spec, weight=weight * bg_params['syn_asymm'])
            nest.Connect(bg_layers['MSN_d2'], bg_layers['MSN_d2'], conn_spec=dict(base_conn_spec, indegree=n_a1), syn_spec=asym_syn)
            nest.Connect(bg_layers['MSN_d2'], bg_layers['MSN_d1'], conn_spec=dict(base_conn_spec, indegree=n_a1), syn_spec=asym_syn)
            nest.Connect(bg_layers['MSN_d1'], bg_layers['MSN_d1'], conn_spec=dict(base_conn_spec, indegree=n_a1), syn_spec=ampa_spec)
            nest.Connect(bg_layers['MSN_d1'], bg_layers['MSN_d2'], conn_spec=dict(base_conn_spec, indegree=n_a2), syn_spec=ampa_spec)

        elif destName == 'MSN':
            nest.Connect(bg_layers[sourceName], bg_layers['MSN_d1'],
                         conn_spec=base_conn_spec, syn_spec=base_syn_spec)
            nest.Connect(bg_layers[sourceName], bg_layers['MSN_d2'],
                         conn_spec=base_conn_spec, syn_spec=base_syn_spec)

        else:
            nest.Connect(bg_layers[sourceName], bg_layers[destName],
                         conn_spec=base_conn_spec, syn_spec=base_syn_spec)

    # Fractional remaining connections: kernel → pairwise_bernoulli
    float_inDegree = inDegree - np.floor(inDegree)
    remaining_connections = np.round(float_inDegree * bg_params['nb' + destName])
    if remaining_connections > 0:
        printv('Adding ' + str(remaining_connections) + ' remaining connections (pairwise_bernoulli)')
        p = 1. / (bg_params['nb' + sourceName] * float(remaining_connections))
        float_conn_spec = {
            'rule': 'pairwise_bernoulli',
            'p': p,
            'mask': {'circular': {'radius': spread}},
            'allow_oversized_mask': True,
            'allow_multapses': True,
        }
        if destName == 'MSN':
            nest.Connect(bg_layers[sourceName], bg_layers['MSN_d1'],
                         conn_spec=float_conn_spec, syn_spec=base_syn_spec)
            nest.Connect(bg_layers[sourceName], bg_layers['MSN_d2'],
                         conn_spec=float_conn_spec, syn_spec=base_syn_spec)
        elif sourceName == 'MSN':
            nest.Connect(bg_layers['MSN_d1'], bg_layers[destName],
                         conn_spec=float_conn_spec, syn_spec=base_syn_spec)
            nest.Connect(bg_layers['MSN_d2'], bg_layers[destName],
                         conn_spec=float_conn_spec, syn_spec=base_syn_spec)
        else:
            nest.Connect(bg_layers[sourceName], bg_layers[destName],
                         conn_spec=float_conn_spec, syn_spec=base_syn_spec)


# -------------------------------------------------------------------------------
# Mirror AMPA connections with NMDA
# NEST 3: SynapseCollection replaces list-of-tuples; use .get() instead of zip(*)
# -------------------------------------------------------------------------------
def mass_mirror_bg(bg_params, nameSrc, nameTgt, source, synapse_label, receptor_type,
                   weight, delay, stochastic_delays=None, verbose=False):
    def printv(text):
        if verbose:
            print(text)

    printv('looking for AMPA connections to mirror with NMDA...\n')
    ampa_conns = nest.GetConnections(source=source, synapse_label=synapse_label)

    if len(ampa_conns) == 0:
        return

    printv('found ' + str(len(ampa_conns)) + ' AMPA connections\n')

    if stochastic_delays is not None and delay > 0:
        printv('Using stochastic delays in mass-mirror')
        delay = np.array(ampa_conns.get('delay')).flatten()

    # NEST 3: extract source/target node IDs from SynapseCollection
    conn_data = ampa_conns.get(['source', 'target'])
    src = nest.NodeCollection(list(conn_data['source']))
    tgt = nest.NodeCollection(list(conn_data['target']))

    if (nameSrc == 'CSN' and nameTgt == 'MSN') or (nameSrc == 'PTN' and nameTgt == 'MSN'):
        if bg_params['plastic_syn']:
            nest.Connect(src, tgt, conn_spec='one_to_one',
                         syn_spec={'synapse_model': 'syn_d1',
                                   'receptor_type': receptor_type,
                                   'weight': weight * bg_params['plast_gain'],
                                   'delay': delay})

            ampa_conns_d2 = nest.GetConnections(source=source, synapse_label=synapse_label + 1000)
            if len(ampa_conns_d2) > 0:
                conn_data_d2 = ampa_conns_d2.get(['source', 'target'])
                src_d2 = nest.NodeCollection(list(conn_data_d2['source']))
                tgt_d2 = nest.NodeCollection(list(conn_data_d2['target']))
                nest.Connect(src_d2, tgt_d2, conn_spec='one_to_one',
                             syn_spec={'synapse_model': 'syn_d2',
                                       'receptor_type': receptor_type,
                                       'weight': weight * bg_params['plast_gain'],
                                       'delay': delay})

            with open('./log/weights_d1.txt', 'a') as f:
                f.write('mirror NMDA d1, weight: ' + str(weight * bg_params['plast_gain']) + ' \n')
                f.write('at mirror connect stdp d1: ' + str(nest.GetDefaults('syn_d1')) + ' \n')
            with open('./log/weights_d2.txt', 'a') as f:
                f.write('mirror NMDA d2, weight: ' + str(weight * bg_params['plast_gain']) + ' \n')
                f.write('at mirror connect stdp d2: ' + str(nest.GetDefaults('syn_d2')) + ' \n')
        else:
            nest.Connect(src, tgt, conn_spec='one_to_one',
                         syn_spec={'synapse_model': 'static_synapse_lbl',
                                   'synapse_label': synapse_label,
                                   'receptor_type': receptor_type,
                                   'weight': weight,
                                   'delay': delay})
    else:
        nest.Connect(src, tgt, conn_spec='one_to_one',
                     syn_spec={'synapse_model': 'static_synapse_lbl',
                               'synapse_label': synapse_label,
                               'receptor_type': receptor_type,
                               'weight': weight,
                               'delay': delay})


def get_frac_bg(bg_params, FractionalOutDegree, nameSrc, nameTgt, cntSrc, cntTgt, useMin=False, verbose=False):
    if not useMin:
        inDegree = get_input_range_bg(bg_params, nameSrc, nameTgt, cntSrc, cntTgt, verbose=verbose)[1] * FractionalOutDegree
    else:
        r = get_input_range_bg(bg_params, nameSrc, nameTgt, cntSrc, cntTgt, verbose=verbose)
        inDegree = (r[1] - r[0]) * FractionalOutDegree + r[0]
    if verbose:
        print('\tConverting fractional outDegree ' + nameSrc + '->' + nameTgt +
              ' from ' + str(FractionalOutDegree) + ' to inDegree: ' + str(round(inDegree, 2)))
    return inDegree


def computeW_bg(bg_params, listRecType, nameSrc, nameTgt, inDegree, gain=1., verbose=False):
    recType = {'AMPA': 1, 'NMDA': 2, 'GABA': 3}
    nu = get_input_range_bg(bg_params, nameSrc, nameTgt,
                            bg_params['count' + nameSrc], bg_params['count' + nameTgt], verbose=verbose)[1]
    if verbose:
        print('\tCompare with the effective chosen inDegree   : ' + str(inDegree))
    LX = bg_params['lx'][nameTgt] * np.sqrt((4. * bg_params['Ri']) / (bg_params['dx'][nameTgt] * bg_params['Rm']))
    attenuation = np.cosh(LX * (1 - bg_params['distcontact'][nameSrc + '->' + nameTgt])) / np.cosh(LX)
    w = {}
    for r in listRecType:
        w[r] = nu / float(inDegree) * attenuation * bg_params['wPSP'][recType[r] - 1] * gain
    return w


def get_input_range_bg(bg_params, nameSrc, nameTgt, cntSrc, cntTgt, verbose=False):
    if nameSrc == 'CSN' or nameSrc == 'PTN':
        nu = bg_params['alpha'][nameSrc + '->' + nameTgt]
        nu0 = 0
        if verbose:
            print('\tMaximal number of distinct input neurons (nu): ' + str(nu))
            print('\tMinimal number of distinct input neurons     : unknown (set to 0)')
    else:
        nu = cntSrc / float(cntTgt) * bg_params['ProjPercent'][nameSrc + '->' + nameTgt] * bg_params['alpha'][nameSrc + '->' + nameTgt]
        nu0 = cntSrc / float(cntTgt) * bg_params['ProjPercent'][nameSrc + '->' + nameTgt]
        if verbose:
            print('\tMaximal number of distinct input neurons (nu): ' + str(nu))
            print('\tMinimal number of distinct input neurons     : ' + str(nu0))
    return [nu0, nu]


def connect_GPi2d_GPi3d(GPi2d, GPi3d):
    # NEST 3: GPi2d and GPi3d are NodeCollections directly
    nest.SetDefaults('static_synapse', {'receptor_type': 0})
    nest.Connect(GPi2d, GPi3d, conn_spec={'rule': 'one_to_one'})
