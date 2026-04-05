#!/usr/bin/env python
# -*- coding: utf-8 -*-
#please set all flags before run the simulation
##

import fetch_params
import ini_all
import nest_routine
import nest
import nest.topology as ntop
import numpy as np
import time
import task
#import gym
import json
import random



def main():
    # 1) reads parameters
    sim_params = fetch_params.read_sim()

    for sim_model in sim_params['sim_model'].keys():
        if sim_params['sim_model'][sim_model]['on']:
            sim_regions=sim_params['sim_model'][sim_model]['regions']
            sim_model_on=sim_model
            print ('simulation model ', sim_model_on, ' will start')
    if sim_regions['S1']:
        ctx_params = fetch_params.read_ctx()
    if sim_regions['M1']:
        ctx_M1_params = fetch_params.read_ctx_M1()
    if sim_regions['TH_S1'] or sim_regions['TH_M1']:
        th_params = fetch_params.read_th()
    if True: #sim_regions['BG']:
        bg_params = fetch_params.read_bg()
    if sim_regions['CB_S1'] or sim_regions['CB_M1']:
        cb_params = fetch_params.read_cb()
    #conn_params = fetch_params.read_conn()

    # 1.5) initialize nest
    nest_routine.initialize_nest(sim_params)

    # 2) instantiates regions
    if sim_regions['S1']:
        ctx_layers = ini_all.instantiate_ctx(ctx_params, sim_params['scalefactor'], sim_params['initial_ignore'])
    if sim_regions['M1']:
        ctx_M1_layers = ini_all.instantiate_ctx_M1(ctx_M1_params, sim_params['scalefactor'],sim_params['initial_ignore'])
    if sim_regions['TH_S1'] or sim_regions['TH_M1']:
        th_layers = ini_all.instantiate_th(th_params, sim_params['scalefactor'],sim_params['initial_ignore'])
    if sim_regions['CB_S1']:
        cb_layers_S1 = ini_all.instantiate_cb('S1', sim_params['scalefactor'], sim_params)
    if sim_regions['CB_M1']:
        cb_layers_M1 = ini_all.instantiate_cb('M1', sim_params['scalefactor'], sim_params)
    if sim_regions['BG']:
        if sim_params['channels']: #if True
            bg_params['channels'] = True #for channels input tasks
        else:
            bg_params['channels'] = False #resting state
        bg_params['circle_center'] = nest_routine.get_channel_centers(sim_params, hex_center=[0, 0],
                                                                ci=sim_params['channels_nb'],
                                                                hex_radius=sim_params['hex_radius'])
        #### Basal ganglia inter-regional connection with S1 and M1 ######
        if sim_regions['S1']:
            if sim_regions['M1']:
                bg_layers,ctx_bg_input = ini_all.instantiate_bg(bg_params, fake_inputs=True,
                                               ctx_inputs={'M1': {'layers': ctx_M1_layers, 'params': ctx_M1_params},
                                                           'S1': {'layers': ctx_layers, 'params': ctx_params},
                                                           'M2': None},
                                                            scalefactor=sim_params['scalefactor'])
            else:
                bg_layers,ctx_bg_input = ini_all.instantiate_bg(bg_params, fake_inputs=True,
                                               ctx_inputs={'M1': None, #{'layers': ctx_M1_layers, 'params': ctx_M1_params},
                                                           'S1': {'layers': ctx_layers, 'params': ctx_params},
                                                           'M2': None},
                                                            scalefactor=sim_params['scalefactor'])
        else:
            if sim_regions['M1']:
                bg_layers,ctx_bg_input = ini_all.instantiate_bg(bg_params, fake_inputs=True,
                                               ctx_inputs={'M1': {'layers': ctx_M1_layers, 'params': ctx_M1_params},
                                                           'S1': None, #{'layers': ctx_layers, 'params': ctx_params},
                                                           'M2': None},
                                                            scalefactor=sim_params['scalefactor'])
            else:
                bg_layers,ctx_bg_input = ini_all.instantiate_bg(bg_params, fake_inputs=True,
                                               ctx_inputs={'M1': None, #{'layers': ctx_M1_layers, 'params': ctx_M1_params},
                                                           'S1': None, #{'layers': ctx_layers, 'params': ctx_params},
                                                           'M2': None},
                                                            scalefactor=sim_params['scalefactor'])

    # 3) interconnect regions
    start_time = time.time()
    if sim_regions['S1'] and sim_regions['CB_S1']:
        _ = nest_routine.connect_region_ctx_cb(ctx_layers['S1_L5B_Pyr'], cb_layers_S1['CB_S1_layer_pons'], 'S1')
    if sim_regions['M1'] and sim_regions['CB_M1']:
        _ = nest_routine.connect_region_ctx_cb(ctx_M1_layers['M1_L5B_PT'], cb_layers_M1['CB_M1_layer_pons'], 'M1')
    if sim_regions['S1'] and sim_regions['TH_S1']:
        _ = nest_routine.connect_region_ctx_th(ctx_layers, th_layers, 'S1')
    if sim_regions['M1'] and sim_regions['TH_M1']:
        _ = nest_routine.connect_region_ctx_th(ctx_M1_layers, th_layers, 'M1')
    if sim_regions['TH_S1'] and sim_regions['S1']:
        _ = nest_routine.connect_region_th_ctx(th_layers, ctx_layers, 'S1')
    if sim_regions['TH_M1'] and sim_regions['M1']:
        _ = nest_routine.connect_region_th_ctx(th_layers, ctx_M1_layers, 'M1')
    if sim_regions['CB_S1'] and sim_regions['S1']:
        _ = nest_routine.connect_region_cb_th(cb_layers_S1, th_layers, 'S1')
    if sim_regions['CB_M1'] and sim_regions['M1']:
        _ = nest_routine.connect_region_cb_th(cb_layers_M1, th_layers, 'M1')
    if sim_regions['BG'] and sim_regions['TH_M1']:
        _ = nest_routine.connect_region_bg_th(bg_layers, th_layers,bg_params)
    with open('./log/' + 'performance.txt', 'a') as file:
        file.write('Interconnect_Regions_Time ' + str(time.time() - start_time) + '\n')

    # 2.5) detectors
    detectors = {}
    if sim_regions['BG']:
        for layer_name in bg_layers.keys():
            detectors[layer_name] = nest_routine.layer_spike_detector(bg_layers[layer_name], layer_name,sim_params['start_time_sp'])
            
            if layer_name=='MSN':
                params={"withgid": True, "withtime": True, "to_file": True, "fbuffer_size":8192}
                params.update({'label': 'MSN_d1', "start": float(sim_params['start_time_sp'])})
                detectors['MSN_d1'] = nest.Create("spike_detector", params= params)
                nest.Connect(pre=nest.GetNodes(bg_layers[layer_name])[0][:int(bg_params['nbMSN']/2)], post=detectors['MSN_d1'])
                
                params.update({'label': 'MSN_d2'})
                detectors['MSN_d2'] = nest.Create("spike_detector", params= params)
                nest.Connect(pre=nest.GetNodes(bg_layers[layer_name])[0][int(bg_params['nbMSN']/2):], post=detectors['MSN_d2'])

    if sim_regions['S1']:
        for layer_name in ctx_layers.keys():
            detectors[layer_name] = nest_routine.layer_spike_detector(ctx_layers[layer_name], layer_name, sim_params['initial_ignore'])
    if sim_regions['M1']:
        for layer_name in ctx_M1_layers.keys():
            detectors[layer_name] = nest_routine.layer_spike_detector(ctx_M1_layers[layer_name], layer_name, sim_params['initial_ignore'])
    if sim_regions['CB_S1']:
        for layer_name in cb_layers_S1.keys():
            detectors[layer_name] = nest_routine.layer_spike_detector(cb_layers_S1[layer_name], layer_name, sim_params['initial_ignore'])
    if sim_regions['CB_M1']:
        for layer_name in cb_layers_M1.keys():
            detectors[layer_name] = nest_routine.layer_spike_detector(cb_layers_M1[layer_name], layer_name, sim_params['initial_ignore'])
    if sim_regions['TH_S1']:
        for layer_name in th_layers['TH_S1_EZ'].keys():
            detectors['TH_S1_EZ' + '_' + layer_name] = nest_routine.layer_spike_detector(th_layers['TH_S1_EZ'][layer_name], 'TH_S1_EZ_'+layer_name, sim_params['initial_ignore'])
        for layer_name in th_layers['TH_S1_IZ'].keys():
            detectors['TH_S1_IZ' + '_' + layer_name] = nest_routine.layer_spike_detector(th_layers['TH_S1_IZ'][layer_name], 'TH_S1_IZ_'+layer_name, sim_params['initial_ignore'])
    if sim_regions['TH_M1']:
        for layer_name in th_layers['TH_M1_EZ'].keys():
            detectors['TH_M1_EZ' + '_' + layer_name] = nest_routine.layer_spike_detector(th_layers['TH_M1_EZ'][layer_name], 'TH_M1_EZ_'+layer_name, sim_params['initial_ignore'])
        for layer_name in th_layers['TH_M1_IZ'].keys():
            detectors['TH_M1_IZ' + '_' + layer_name] = nest_routine.layer_spike_detector(th_layers['TH_M1_IZ'][layer_name], 'TH_M1_IZ_'+layer_name, sim_params['initial_ignore'])
    print (sim_model_on)

    if sim_model_on=='resting_state':
        #nest_routine.get_connections_to_file('MSN','GPe',bg_layers['MSN'],bg_layers['GPe'])
        #nest_routine.get_connections_to_file('MSN','GPi',bg_layers['MSN'],bg_layers['GPi'])

        simulation_time = sim_params['simDuration']+sim_params['initial_ignore']
        print('Simulation Started:')
        start_time = time.time()
        nest.Simulate(simulation_time)
        with open('./log/' + 'performance.txt', 'a') as file:
            file.write('Simulation_Elapse_Time ' + str(time.time() - start_time) + '\n')
        
        mean_fr = {}
        At = {}
        if sim_regions['BG']:
          for layer_name in list(bg_layers.keys())+['MSN_d1','MSN_d2']:
            if layer_name == 'MSN_d1' or layer_name == 'MSN_d2':
                rate = nest_routine.average_fr(detectors[layer_name], sim_params['simDuration']-sim_params['start_time_sp'],int(nest_routine.count_layer(bg_layers[layer_name[:3]])/2))
                act = nest_routine.instantaneous_fr(detectors[layer_name], sim_params['start_time_sp'], sim_params['simDuration'],int(nest_routine.count_layer(bg_layers[layer_name[:3]])/2),float(sim_params['dt']))
            else:
                rate = nest_routine.average_fr(detectors[layer_name], sim_params['simDuration']-sim_params['start_time_sp'],nest_routine.count_layer(bg_layers[layer_name]))
                act = nest_routine.instantaneous_fr(detectors[layer_name], sim_params['start_time_sp'], sim_params['simDuration'],nest_routine.count_layer(bg_layers[layer_name]),float(sim_params['dt']))
            mean_fr[layer_name]=rate
            At[layer_name]=list(act)
        with open('./log/mean_fr.json','w') as my_file:
          json.dump(mean_fr, my_file)
        with open('./log/At.json','w') as my_At_file:
          json.dump(At, my_At_file)

        print('end resting')

    elif sim_model_on == 'action_selection':
        ###########################################################################################################################
        ################### Simulated Input on PTN and CSN ##########################################
        ###########################################################################################################################
        sample_input_size = 200 
        ##### get columns data , create PG and connect ############
        CSN_columns_gids = nest_routine.get_columns_data('CSN', bg_params['circle_center'],sim_params['channels_radius']) 
        CSN_psg, CSN_syn = nest_routine.create_psg_channels(syn_weight=1., syn_delay=2.,channels_nb=sim_params['channels_nb']) #0.3 #tested: 27.5 (looks caothic), 7.5, 3., 1., 0.5, 0.1, 5.#original: 11.5
        #select a subset of 500 neurons for each channel
        CSN_columns_gids_subset = []
        for k1 in CSN_columns_gids:
            CSN_columns_gids_subset.append(random.sample(k1,sample_input_size))
        #_ = nest_routine.connect_psg_to_channels(CSN_columns_gids, CSN_psg, CSN_syn)
        _ = nest_routine.connect_psg_to_channels(CSN_columns_gids_subset, CSN_psg, CSN_syn)

        PTN_columns_gids = nest_routine.get_columns_data('PTN', bg_params['circle_center'],sim_params['channels_radius']) 
        PTN_psg, PTN_syn = nest_routine.create_psg_channels(syn_weight=1., syn_delay=2.,channels_nb=sim_params['channels_nb']) #0.3 #tested: 27.5 (looks caothic), 7.5, 3., 1., 0.5, 0.1, 5.#original: 11.5
        #select a subset of 500 neurons for each channel
        PTN_columns_gids_subset = []
        for k2 in PTN_columns_gids:
            PTN_columns_gids_subset.append(random.sample(k2,sample_input_size))
        #_ = nest_routine.connect_psg_to_channels(PTN_columns_gids, PTN_psg, PTN_syn)
        _ = nest_routine.connect_psg_to_channels(PTN_columns_gids_subset, PTN_psg, PTN_syn)

        MSN_columns_gids = nest_routine.get_columns_data('MSN', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN
        MSN_d1_columns_gids = nest_routine.get_columns_data('MSN_d1', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN
        MSN_d2_columns_gids = nest_routine.get_columns_data('MSN_d2', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN

        GPi_columns_gids = nest_routine.get_columns_data('GPi', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of GPi

        channel_1_idx = 1
        channel_2_idx = 2 

        nest.Simulate(1000.)
        stim_dt = 500.
        stim_break = 500.
        episode = 0
        MSN_fr, MSN_d1_fr, MSN_d2_fr, GPi_fr, PTN_fr, CSN_fr = [],[],[],[], [], []
        PTN_CSN_ratio = 1.8 #2.

        ##### download wiring for CTX - MSN ######
        gid_pos_MSN_d1 = np.loadtxt('./log/'+'MSN_d1'+'.txt')
        gid_pos_MSN_d2 = np.loadtxt('./log/'+'MSN_d2'+'.txt')
        for ni,i in enumerate(CSN_columns_gids_subset): #where those 200 neurons are connected ?
            np.savetxt('./log/'+'CSN_c'+str(ni)+'_to_'+'MSN_d1'+'.txt',nest.GetStatus(nest.GetConnections(source=[gid[0] for gid in i],target=list(gid_pos_MSN_d1[:,0].astype(int))),keys={'source','target','weight'}))
            np.savetxt('./log/'+'CSN_c'+str(ni)+'_to_'+'MSN_d2'+'.txt',nest.GetStatus(nest.GetConnections(source=[gid[0] for gid in i],target=list(gid_pos_MSN_d2[:,0].astype(int))),keys={'source','target','weight'}))
        for ni,i in enumerate(PTN_columns_gids_subset): #where those 200 neurons are connected ?
            np.savetxt('./log/'+'PTN_c'+str(ni)+'_to_'+'MSN_d1'+'.txt',nest.GetStatus(nest.GetConnections(source=[gid[0] for gid in i],target=list(gid_pos_MSN_d1[:,0].astype(int))),keys={'source','target','weight'}))
            np.savetxt('./log/'+'PTN_c'+str(ni)+'_to_'+'MSN_d2'+'.txt',nest.GetStatus(nest.GetConnections(source=[gid[0] for gid in i],target=list(gid_pos_MSN_d2[:,0].astype(int))),keys={'source','target','weight'}))
        ##########################################

        for i in np.arange(11): #11): #channel 1
            for j in np.arange(11): #11):  #channel 2

                Poisson_ch1 = i*1.7 #i*1.2 #i*0.9 #i*0.8 #2.15 #0.2 
                Poisson_ch2 = j*1.7 #j*1.2 #j*0.9 #j*0.8 #2.15 #0.2
                
                stim_start = 1000. + (stim_dt+stim_break)*episode
                stim_stop =  stim_start + stim_dt

                #stimulus on channel 1
                _ = nest_routine.apply_fixed_stimulus(time_start=stim_start, time_stop=stim_stop, psg=CSN_psg[channel_1_idx], stimulus=Poisson_ch1) #M1
                _ = nest_routine.apply_fixed_stimulus(time_start=stim_start, time_stop=stim_stop, psg=PTN_psg[channel_1_idx], stimulus=Poisson_ch1*PTN_CSN_ratio)
                #stimulus on channel 2
                _ = nest_routine.apply_fixed_stimulus(time_start=stim_start, time_stop=stim_stop, psg=CSN_psg[channel_2_idx], stimulus=Poisson_ch2) #M1
                _ = nest_routine.apply_fixed_stimulus(time_start=stim_start, time_stop=stim_stop, psg=PTN_psg[channel_2_idx], stimulus=Poisson_ch2*PTN_CSN_ratio)

                nest.Simulate(stim_dt + stim_break)

                print('Simulating episode ',episode,'  stimulation from:',stim_start,'  to:',stim_stop)

                ####### For MSN (get mean firing rate )####################################
                sp_MSN = nest.GetStatus(detectors['MSN'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_MSN = np.where(np.logical_and(sp_MSN['times']>=stim_start, sp_MSN['times']<=stim_stop))
                print('MSN')
                MSN_fr.append(nest_routine.get_targets_mean_rates(MSN_columns_gids,sp_MSN,idx_dt_MSN,stim_start,stim_stop))

                sp_MSN_d1 = nest.GetStatus(detectors['MSN_d1'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_MSN_d1 = np.where(np.logical_and(sp_MSN_d1['times']>=stim_start, sp_MSN_d1['times']<=stim_stop))
                print('MSN_d1')
                MSN_d1_fr.append(nest_routine.get_targets_mean_rates(MSN_d1_columns_gids,sp_MSN_d1,idx_dt_MSN_d1,stim_start,stim_stop))

                sp_MSN_d2 = nest.GetStatus(detectors['MSN_d2'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_MSN_d2 = np.where(np.logical_and(sp_MSN_d2['times']>=stim_start, sp_MSN_d2['times']<=stim_stop))
                print('MSN_d2')
                MSN_d2_fr.append(nest_routine.get_targets_mean_rates(MSN_d2_columns_gids,sp_MSN_d2,idx_dt_MSN_d2,stim_start,stim_stop))

                ####### For GPi (get mean firing rate ) ####################################
                sp_GPi = nest.GetStatus(detectors['GPi'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_GPi = np.where(np.logical_and(sp_GPi['times']>=stim_start, sp_GPi['times']<=stim_stop))  #M1
                print('GPi')
                GPi_fr.append(nest_routine.get_targets_mean_rates(GPi_columns_gids,sp_GPi,idx_dt_GPi,stim_start,stim_stop)) #M1
                ####### For MSN (get mean firing rate )####################################
                sp_PTN = nest.GetStatus(detectors['PTN'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_PTN = np.where(np.logical_and(sp_PTN['times']>=stim_start, sp_PTN['times']<=stim_stop))
                print('PTN')
                #PTN_fr.append(nest_routine.get_targets_mean_rates(PTN_columns_gids,sp_PTN,idx_dt_PTN,stim_start,stim_stop))
                PTN_fr.append(nest_routine.get_targets_mean_rates(PTN_columns_gids_subset,sp_PTN,idx_dt_PTN,stim_start,stim_stop))
                ####### For MSN (get mean firing rate )####################################
                sp_CSN = nest.GetStatus(detectors['CSN'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_CSN = np.where(np.logical_and(sp_CSN['times']>=stim_start, sp_CSN['times']<=stim_stop))
                print('CSN')
                #CSN_fr.append(nest_routine.get_targets_mean_rates(CSN_columns_gids,sp_CSN,idx_dt_CSN,stim_start,stim_stop))
                CSN_fr.append(nest_routine.get_targets_mean_rates(CSN_columns_gids_subset,sp_CSN,idx_dt_CSN,stim_start,stim_stop))

                with open('./log/' + 'GPi_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in GPi_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'MSN_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in MSN_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'MSN_d1_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in MSN_d1_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'MSN_d2_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in MSN_d2_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'PTN_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in PTN_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'CSN_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in CSN_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')


                episode+=1
        
        #save data for analysis#
        #np.savetxt('./log/GPi_fr.txt',GPi_fr)
        #np.savetxt('./log/MSN_fr.txt',MSN_fr)
        #np.savetxt('./log/MSN_d1_fr.txt',MSN_d1_fr)
        #np.savetxt('./log/MSN_d2_fr.txt',MSN_d2_fr)
        #np.savetxt('./log/PTN_fr.txt',PTN_fr)
        #np.savetxt('./log/CSN_fr.txt',CSN_fr)

        print('finished')
   
    elif sim_model_on == 'plasticity':
        if bg_params['plastic_syn']:

            DA_spikes_d1  =nest.Create('poisson_generator')
            my_parrot_d1 = nest.Create('parrot_neuron')
            nest.Connect(DA_spikes_d1,my_parrot_d1)
            #nest.Connect(dopa_neurons,my_parrot)
            nest.Connect(my_parrot_d1,bg_params['vt_d1']) #syn_weight=1.

            DA_spikes_d2  =nest.Create('poisson_generator')
            my_parrot_d2 = nest.Create('parrot_neuron')
            nest.Connect(DA_spikes_d2,my_parrot_d2)
            #nest.Connect(dopa_neurons,my_parrot)
            nest.Connect(my_parrot_d2,bg_params['vt_d2']) #syn_weight=1.

            ##### get columns data , create PG and connect ############
            sample_input_size = 300 #200 

            CSN_columns_gids = nest_routine.get_columns_data('CSN', bg_params['circle_center'],sim_params['channels_radius']) 
            CSN_psg, CSN_syn = nest_routine.create_psg_channels(syn_weight=1., syn_delay=2.,channels_nb=sim_params['channels_nb']) #0.3 #tested: 27.5 (looks caothic), 7.5, 3., 1., 0.5, 0.1, 5.#original: 11.5
            #select a subset of 500 neurons for each channel
            CSN_columns_gids_subset = []
            for k1 in CSN_columns_gids:
                CSN_columns_gids_subset.append(random.sample(k1,sample_input_size))
            #_ = nest_routine.connect_psg_to_channels(CSN_columns_gids, CSN_psg, CSN_syn)
            _ = nest_routine.connect_psg_to_channels(CSN_columns_gids_subset, CSN_psg, CSN_syn)
            
            PTN_columns_gids = nest_routine.get_columns_data('PTN', bg_params['circle_center'],sim_params['channels_radius']) 
            PTN_psg, PTN_syn = nest_routine.create_psg_channels(syn_weight=1., syn_delay=2.,channels_nb=sim_params['channels_nb']) #0.3 #tested: 27.5 (looks caothic), 7.5, 3., 1., 0.5, 0.1, 5.#original: 11.5
            #select a subset of 500 neurons for each channel
            PTN_columns_gids_subset = []
            for k2 in PTN_columns_gids:
                PTN_columns_gids_subset.append(random.sample(k2,sample_input_size))
            #_ = nest_routine.connect_psg_to_channels(PTN_columns_gids, PTN_psg, PTN_syn)
            _ = nest_routine.connect_psg_to_channels(PTN_columns_gids_subset, PTN_psg, PTN_syn)
            
            MSN_columns_gids = nest_routine.get_columns_data('MSN', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN
            MSN_d1_columns_gids = nest_routine.get_columns_data('MSN_d1', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN
            MSN_d2_columns_gids = nest_routine.get_columns_data('MSN_d2', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN

            GPi_columns_gids = nest_routine.get_columns_data('GPi', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of GPi
            
            #CMPf_columns_gids = nest_routine.get_columns_data('CMPf', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of CMPf
            FSI_columns_gids = nest_routine.get_columns_data('FSI', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of FSI
            GPe_columns_gids = nest_routine.get_columns_data('GPe', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of GPe
            STN_columns_gids = nest_routine.get_columns_data('STN', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of STN

            channel_1_idx = 1
            channel_2_idx = 2 

            global modulated
            modulated = False


            nest.Simulate(1000.)
            stim_dt = 1000.#500.
            stim_break = 500.#1000. #time after stimulation offset
            dopa_break = 1000.#500. #time after dopamine release 
            dopa_dt = 300.#10.
            MSN_fr, GPi_fr, PTN_fr, CSN_fr, MSN_d1_fr, MSN_d2_fr, CMPf_fr, FSI_fr, GPe_fr, STN_fr = [],[],[],[],[],[],[],[],[],[]
            CS_plus, CS_minus = [],[]
            
            Poisson_ch1 = 5.*1.7 #5*1.7 #4*1.7 #6*1.7 #6*0.8 #7*0.8 #0.2 
            Poisson_ch2 = 2.*1.7 #2*1.7 #4*1.7 #5*0.8 #5*0.8 #0.2
            dopa_rate = 20.0 #23.0 #25.0 #30.0 #25.0 #20.0 #25. #30. #10. #80. #10. #80.

            PTN_CSN_ratio = 1.8 #2.
            hits = 0

            #cond_from = 600 #700 #600 #400 #60 #80
            #total_episodes = 600 #700 #600 #400 #150

            current_task = [] #1:CS+, -1:CS-
            con_ep = 100 #150 #70 #150
            dis_ep= 100 #150  #70 #150
            gen_ep= 30   #50
            cond_from = con_ep + dis_ep + gen_ep
            total_episodes = cond_from
            dc = 0.0042 #0.004 #0.0038 #0.0035 #0.004 #0.0045 #0.0055 #0.006 #0.007 #0.012 #0.013   

            for episode in np.arange(total_episodes):
                if episode > cond_from: #lets make 80*2 CS- , and 40*2 CS+
                    if random.random() > 1.0: #0.97: #0.7: ## this is CS+
                        CS_plus.append(1)
                        CS_minus.append(0)
                        Poisson_ch1 = 5*1.7 #4*1.7 #7*0.8 #0.2 
                        Poisson_ch2 = 2*1.7 #5*0.8 #0.2
                        current_task.append(1)
                    else:  ##this is CS-
                        CS_minus.append(1)
                        CS_plus.append(0)
                        Poisson_ch1 = 2*1.7 #7*0.8 #0.2 
                        Poisson_ch2 = 5*1.7 #4*1.7 #5*0.8 #0.2
                        current_task.append(-1)
                        #STDPs for DA dips
                        nest.SetDefaults('syn_d1',{'vt':bg_params['vt_d2'][0],
                             'A_plus':-1.0*dc/40.0,'A_minus':dc/40.0}) 
                        nest.SetDefaults('syn_d2',{'vt':bg_params['vt_d2'][0],
                             'A_plus':dc,'A_minus':-1.0*dc})

                else:
                    CS_plus.append(1)
                    CS_minus.append(0)
                    current_task.append(1)
                    #STDPs for DA bursts
                    nest.SetDefaults('syn_d1',{'vt':bg_params['vt_d1'][0],
                             'A_plus':dc,'A_minus':(dc/4.0)}) 
                    nest.SetDefaults('syn_d2',{'vt':bg_params['vt_d1'][0],
                             'A_plus':dc/3.0,'A_minus':dc*3.0/4.0})

                #if episode==96: #total_episodes-2: #one before the last
                #    Poisson_ch1 = 2*1.7 #7*0.8 #0.2 
                #    Poisson_ch2 = 5*1.7 #4*1.7 #5*0.8 #0.2
                #if episode==100: #total_episodes-1: #the last one
                #    Poisson_ch1 = 5*1.7 #4*1.7 #7*0.8 #0.2 
                #    Poisson_ch2 = 2*1.7 #5*0.8 #0.2
                
                if episode > ( con_ep + gen_ep):
                    with open('./log/' + 'do_tests.txt', 'r') as file:
                        for lines in file:
                            do_tests = lines
                    print('do_tests: ',do_tests,'length: ',len(do_tests))
                    if do_tests[:-1]=='T2':
                        Poisson_ch1 = 2*1.7 #7*0.8 #0.2
                        Poisson_ch2 = 5*1.7 #4*1.7 #5*0.8 #0.2
                    if do_tests[:-1]=='T1': #total_episodes-1: #the last one
                        Poisson_ch1 = 5*1.7 #4*1.7 #7*0.8 #0.2
                        Poisson_ch2 = 2*1.7 #5*0.8 #0.2
                else:
                    if episode < (con_ep + 1):
                        do_tests = 'F '
                    else:
                        do_tests = 'T '
                        Poisson_ch1 = 0.0
                        Poisson_ch2 = 0.0
                    ##############
                    if do_tests[:-1]=='T' and (episode==(con_ep + 5) or episode==(con_ep + 15)):
                        do_tests = 'T1 '
                        Poisson_ch1 = 5*1.7 #4*1.7 #7*0.8 #0.2
                        Poisson_ch2 = 2*1.7 #5*0.8 #0.2

                    if do_tests[:-1]=='T' and (episode==(con_ep + 10)  or episode==(con_ep + 20)):
                        do_tests= 'T2 '
                        Poisson_ch1 = 2*1.7 #7*0.8 #0.2
                        Poisson_ch2 = 5*1.7 #4*1.7 #5*0.8 #0.2
                    ##############

                stim_start = 1000. + (stim_dt + stim_break + dopa_dt + dopa_break)*episode
                stim_stop =  stim_start + stim_dt
                dopa_start = stim_stop +  stim_break
                dopa_end = dopa_start + dopa_dt

                #stimulus on channel 1
                _ = nest_routine.apply_fixed_stimulus(time_start=stim_start, time_stop=stim_stop, psg=CSN_psg[channel_1_idx], stimulus=Poisson_ch1) #M1
                _ = nest_routine.apply_fixed_stimulus(time_start=stim_start, time_stop=stim_stop, psg=PTN_psg[channel_1_idx], stimulus=Poisson_ch1*PTN_CSN_ratio)
                #stimulus on channel 2
                _ = nest_routine.apply_fixed_stimulus(time_start=stim_start, time_stop=stim_stop, psg=CSN_psg[channel_2_idx], stimulus=Poisson_ch2) #M1
                _ = nest_routine.apply_fixed_stimulus(time_start=stim_start, time_stop=stim_stop, psg=PTN_psg[channel_2_idx], stimulus=Poisson_ch2*PTN_CSN_ratio)

                #Simulation
                nest.Simulate(stim_dt + stim_break)

                print('Stimulation from:',stim_start,'  to:',stim_stop,'\n')
                ####### For MSN (get mean firing rate )####################################
                sp_MSN = nest.GetStatus(detectors['MSN'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_MSN = np.where(np.logical_and(sp_MSN['times']>=stim_start, sp_MSN['times']<=stim_stop))
                print('MSN')
                MSN_fr.append(nest_routine.get_targets_mean_rates(MSN_columns_gids,sp_MSN,idx_dt_MSN,stim_start,stim_stop))
                
                sp_MSN_d1 = nest.GetStatus(detectors['MSN_d1'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_MSN_d1 = np.where(np.logical_and(sp_MSN_d1['times']>=stim_start, sp_MSN_d1['times']<=stim_stop))
                print('MSN_d1')
                MSN_d1_fr.append(nest_routine.get_targets_mean_rates(MSN_d1_columns_gids,sp_MSN_d1,idx_dt_MSN_d1,stim_start,stim_stop))
                
                sp_MSN_d2 = nest.GetStatus(detectors['MSN_d2'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_MSN_d2 = np.where(np.logical_and(sp_MSN_d2['times']>=stim_start, sp_MSN_d2['times']<=stim_stop))
                print('MSN_d2')
                MSN_d2_fr.append(nest_routine.get_targets_mean_rates(MSN_d2_columns_gids,sp_MSN_d2,idx_dt_MSN_d2,stim_start,stim_stop))

                ####### For GPi (get mean firing rate ) ####################################
                sp_GPi = nest.GetStatus(detectors['GPi'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_GPi = np.where(np.logical_and(sp_GPi['times']>=stim_start, sp_GPi['times']<=stim_stop))  #M1
                print('GPi')
                GPi_fr.append(nest_routine.get_targets_mean_rates(GPi_columns_gids,sp_GPi,idx_dt_GPi,stim_start,stim_stop)) #M1
                ####### For MSN (get mean firing rate )####################################
                #sp_PTN = nest.GetStatus(detectors['PTN'],'events')[0] # will get data for MSN_columns_gids
                #idx_dt_PTN = np.where(np.logical_and(sp_PTN['times']>=stim_start, sp_PTN['times']<=stim_stop))
                #print('PTN')
                #PTN_fr.append(nest_routine.get_targets_mean_rates(PTN_columns_gids_subset,sp_PTN,idx_dt_PTN,stim_start,stim_stop))
                ####### For MSN (get mean firing rate )####################################
                #sp_CSN = nest.GetStatus(detectors['CSN'],'events')[0] # will get data for MSN_columns_gids
                #idx_dt_CSN = np.where(np.logical_and(sp_CSN['times']>=stim_start, sp_CSN['times']<=stim_stop))
                #print('CSN')
                #CSN_fr.append(nest_routine.get_targets_mean_rates(CSN_columns_gids_subset,sp_CSN,idx_dt_CSN,stim_start,stim_stop))
                ###### get the rest of the firing rates: ###############
                #sp_CMPf = nest.GetStatus(detectors['CMPf'],'events')[0] # will get data for MSN_columns_gids
                #idx_dt_CMPf = np.where(np.logical_and(sp_CMPf['times']>=stim_start, sp_CMPf['times']<=stim_stop))
                #print('CMPf')
                #CMPf_fr.append(nest_routine.get_targets_mean_rates(CMPf_columns_gids,sp_CMPf,idx_dt_CMPf,stim_start,stim_stop))

                sp_FSI = nest.GetStatus(detectors['FSI'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_FSI = np.where(np.logical_and(sp_FSI['times']>=stim_start, sp_FSI['times']<=stim_stop))
                print('FSI')
                FSI_fr.append(nest_routine.get_targets_mean_rates(FSI_columns_gids,sp_FSI,idx_dt_FSI,stim_start,stim_stop))

                sp_GPe = nest.GetStatus(detectors['GPe'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_GPe = np.where(np.logical_and(sp_GPe['times']>=stim_start, sp_GPe['times']<=stim_stop))
                print('GPe')
                GPe_fr.append(nest_routine.get_targets_mean_rates(GPe_columns_gids,sp_GPe,idx_dt_GPe,stim_start,stim_stop))

                sp_STN = nest.GetStatus(detectors['STN'],'events')[0] # will get data for MSN_columns_gids
                idx_dt_STN = np.where(np.logical_and(sp_STN['times']>=stim_start, sp_STN['times']<=stim_stop))
                print('STN')
                STN_fr.append(nest_routine.get_targets_mean_rates(STN_columns_gids,sp_STN,idx_dt_STN,stim_start,stim_stop))


                ## Task evaluation ##
                '''
                if GPi_fr[-1][channel_1_idx] < GPi_fr[-1][channel_2_idx]: #if channel 1 is winner
                    #activation of VTA+SNc
                    if current_task[-1] == 1: #CS+ task
                        _ = nest.SetStatus(DA_spikes_d1,{'rate':dopa_rate,'start':dopa_start,'stop':dopa_end}) #rate: 150., 20 100., 200., 300.
                        print('episode: ',episode,', Dopamine burst from vt_d1:',dopa_start,'  to:',dopa_end,'\n')
                else:
                    if episode > cond_from: #lets make 80*2 CS- , and 40*2 CS+
                        if current_task[-1] == -1: #CS- task
                            _ = nest.SetStatus(DA_spikes_d2,{'rate':dopa_rate,'start':dopa_start,'stop':dopa_end+100.}) #+ 100 ms for the case of dip
                            print('episode: ',episode,', Dopamine Release from vt_d2:',dopa_start,'  to:',dopa_end+100.,'\n')
                '''     

                if episode < total_episodes-2 and do_tests[:-1]=='F':  #95:

                    if current_task[-1] == 1: #CS+ task
                        if GPi_fr[-1][channel_1_idx] < GPi_fr[-1][channel_2_idx]: #if channel 1 is winner
                            _ = nest.SetStatus(DA_spikes_d1,{'rate':dopa_rate,'start':dopa_start,'stop':dopa_end}) #rate: 150., 20 100., 200., 300.
                            print('episode: ',episode,' ,task: ',current_task[-1],', Dopamine burst from vt_d1:',dopa_start,'  to:',dopa_end,'\n')
                
                    if current_task[-1] == -1 and MSN_d2_fr[-1][2]<70.0: #0.7: #CS- task and episode > cond_from: #lets make 80*2 CS- , and 40*2 CS+
                        if GPi_fr[-1][channel_1_idx] < GPi_fr[-1][channel_2_idx]: #if channel 1 is winner
                            _ = nest.SetStatus(DA_spikes_d2,{'rate':dopa_rate,'start':dopa_start,'stop':dopa_end+100.}) #+ 100 ms for the case of dip
                            print('episode: ',episode,' ,task: ',current_task[-1],', Dopamine Release from vt_d2:',dopa_start,'  to:',dopa_end+100.,'\n')
                            
                            #modulation of D2 neurons synapses.
                            if not(modulated): #then modulate.
                                my_gain = bg_params['MSND2_mod'] #2.5

                                my_target = nest.GetNodes(bg_layers['MSN'])[0][:]#[:int(bg_params['nbMSN']/2)] #'MSN_d1'
                                my_source = nest.GetNodes(bg_layers['MSN'])[0][int(bg_params['nbMSN']/2):] #'MSN_d2'
                                conn = nest.GetConnections(source=my_source, target=my_target)
                                old_weight = nest.GetStatus(conn)[0]['weight']
                                #print(nest.GetStatus(conn))
                                newweight = old_weight*my_gain #2.5 #2.0 #4.0 #2.0
                                nest.SetStatus(conn, {'weight': newweight})

                                my_target_1 = nest.GetNodes(bg_layers['GPe'])[0]
                                conn_1 = nest.GetConnections(source=my_source, target=my_target_1)
                                old_weight_1 = nest.GetStatus(conn_1)[0]['weight']
                                #print(nest.GetStatus(conn))
                                newweight_1 = old_weight_1*my_gain #2.0 #4.0 #2.0
                                nest.SetStatus(conn_1, {'weight': newweight_1})

                                my_target_2 = nest.GetNodes(bg_layers['GPi'])[0]
                                conn_2 = nest.GetConnections(source=my_source, target=my_target_2)
                                old_weight_2 = nest.GetStatus(conn_2)[0]['weight']
                                #print(nest.GetStatus(conn))
                                newweight_2 = old_weight_2*my_gain #2.0 #4.0 #2.0
                                nest.SetStatus(conn_2, {'weight': newweight_2})

                                modulated = True
                                #print(nest.GetStatus(conn))
                        else:
                            print('episode: ',episode,' ,task: ',current_task[-1],'no dopamine release, selected channel 2')
                else:
                    print('testing GPi response without Dop burst or dip')


                #if MSN_d1_fr[-1][1]>4.0: #1.0: #for channel 1 
                if episode == con_ep: #150 : # we will change to CS- after 200 episodes. 
                    cond_from = episode
                    #stim_break = 100.0


                nest.Simulate(dopa_dt + dopa_break)

                #Record results#
                with open('./log/' + 'GPi_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in GPi_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'MSN_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in MSN_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'MSN_d1_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in MSN_d1_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'MSN_d2_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in MSN_d2_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                #with open('./log/' + 'PTN_fr.txt', 'a') as file:
                #    r_val_chain = ''
                #    for r_val in PTN_fr[-1]:
                #        r_val_chain = r_val_chain + str(r_val) + ' '
                #    file.write(r_val_chain[:-1]+ '\n')
                #with open('./log/' + 'CSN_fr.txt', 'a') as file:
                #    r_val_chain = ''
                #    for r_val in CSN_fr[-1]:
                #        r_val_chain = r_val_chain + str(r_val) + ' '
                #    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'FSI_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in FSI_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'GPe_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in GPe_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')
                #with open('./log/' + 'CMPf_fr.txt', 'a') as file:
                #    r_val_chain = ''
                #    for r_val in CMPf_fr[-1]:
                #        r_val_chain = r_val_chain + str(r_val) + ' '
                #    file.write(r_val_chain[:-1]+ '\n')
                with open('./log/' + 'STN_fr.txt', 'a') as file:
                    r_val_chain = ''
                    for r_val in STN_fr[-1]:
                        r_val_chain = r_val_chain + str(r_val) + ' '
                    file.write(r_val_chain[:-1]+ '\n')

                with open('./log/' + 'CS_plus.txt', 'a') as file:
                    file.write(str(CS_plus[-1])+ '\n')

                with open('./log/' + 'CS_minus.txt', 'a') as file:
                    file.write(str(CS_minus[-1])+ '\n')
                
                with open('./log/' + 'current_task.txt', 'a') as file:
                    file.write(str(current_task[-1])+ '\n')


            #save data for analysis#
            #np.savetxt('./log/GPi_fr.txt',GPi_fr)
            #np.savetxt('./log/MSN_fr.txt',MSN_fr)
            #np.savetxt('./log/MSN_d1_fr.txt',MSN_d1_fr)
            #np.savetxt('./log/MSN_d2_fr.txt',MSN_d2_fr)
            #np.savetxt('./log/PTN_fr.txt',PTN_fr)
            #np.savetxt('./log/CSN_fr.txt',CSN_fr)
            #np.savetxt('./log/CS_plus.txt',CS_plus)
            #np.savetxt('./log/CS_minus.txt',CS_minus)

            print('finished')

        else:
            print('error: plastic_syn flag is deactivated')
  
    elif sim_model_on == 'reinf_learning':
        '''
        ##################################################################################################################
        ################### simulated Input on S1 L5B and L5A Pyr from hand position by Poisson generator ################
        ##################################################################################################################
        ## get target GIDS 
        S1_L5A_columns_gids = nest_routine.get_columns_data('S1_L5A_Pyr', bg_params['circle_center'],sim_params['channels_radius'])
        S1_L5B_columns_gids = nest_routine.get_columns_data('S1_L5B_Pyr', bg_params['circle_center'],sim_params['channels_radius'])
        print('S1 L5B columns gids: ',S1_L5B_columns_gids)
        ## create stimulus based on distance. They are diccionaries to be used in function set_stimulus ####
        stimulus_L5A = {'scaling_stim':30*17.7, 'bias':2.,'mu':0.,'sig': 0.25} 
        stimulus_L5B = {'scaling_stim':30*31.3, 'bias':15.,'mu':0.,'sig':0.25} 
        ## create PG and synapses
        S1_L5A_psg, S1_L5A_syn = nest_routine.create_psg_channels(syn_weight=1.5, syn_delay=1.5,channels_nb=sim_params['channels_nb'])  #tested: 4. #original: 2.
        S1_L5B_psg, S1_L5B_syn = nest_routine.create_psg_channels(syn_weight=1.7, syn_delay=1.5,channels_nb=sim_params['channels_nb']) #tested:4.2 #original: 2.2
        ## connect PG to GIDS
        _ = nest_routine.connect_psg_to_channels(S1_L5A_columns_gids, S1_L5A_psg, S1_L5A_syn)
        _ = nest_routine.connect_psg_to_channels(S1_L5B_columns_gids, S1_L5B_psg, S1_L5B_syn)
        '''


        ##################################### get targets Ids ##################################################################
        ## get target GIDS
        MSN_columns_gids = nest_routine.get_columns_data('MSN', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN
        GPi_columns_gids = nest_routine.get_columns_data('GPi', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of GPi
        TH_M1_IZ_TC_columns_gids = nest_routine.get_columns_data('TH_M1_IZ_thalamic_nucleus_TC', bg_params['circle_center'],sim_params['channels_radius'])#0.05)#sim_params['channels_radius'])  #this gids are used for calculate mean rate of GPi
    
        
        ###########################################################################################################################
        ################### Simulated Input on MSN from M2 by Poisson generator M2 = CSN ##########################################
        ###########################################################################################################################
        ##### get columns data , create PG and connect ############
        CSN_columns_gids = nest_routine.get_columns_data('CSN', bg_params['circle_center'],sim_params['channels_radius']) 
        stimulus_CSN = {'scaling_stim':3.3, 'bias':0.,'mu':0.,'sig':0.25} #scaling_stim':3.3, 0.5 #2.3 #nest_routine.create_stimulus(scaling_stim=3.45, bias=0.)#3.5 looks ok, maybe the best #3.3 ok, but a little low, 3.5 produces spikes in near-by channels, looks good-> 4. #1. with weight 27.5 shows no response at MSN #10. gave around 30 Hz into the MSN
        CSN_psg, CSN_syn = nest_routine.create_psg_channels(syn_weight=0.3, syn_delay=1.5,channels_nb=sim_params['channels_nb']) #0.3 #tested: 27.5 (looks caothic), 7.5, 3., 1., 0.5, 0.1, 5.#original: 11.5
        _ = nest_routine.connect_psg_to_channels(CSN_columns_gids, CSN_psg, CSN_syn)


        ######################################################################################################################
        ################################################## PTN "Hand position" as S1 ################################################################
        ##### create a grid of references points in the grid #########
        grid_centers = nest_routine.get_sheet_grid(10) #5 #5+1 x 5+1 (36 elements distribuited uniformly on the 2D sheet)
        ##### get columns data , create PG and connect ############
        PTN_columns_gids = nest_routine.get_columns_data('PTN', grid_centers,0.05)#sim_params['channels_radius']) #ok
        #stimulus_PTN = {'scaling_stim':4*17.7, 'bias':2.,'mu':0.,'sig': 0.25} 
        PTN_psg, PTN_syn = nest_routine.create_psg_channels(syn_weight=0.0001, syn_delay=1.5,channels_nb=len(grid_centers)) #tested:4.2 #original: 2.2
        _ = nest_routine.connect_psg_to_channels(PTN_columns_gids, PTN_psg, PTN_syn)
        rate_Poisson = 8.#2.#10.#50.#100.#200.#500.


        ################### Devices for dopamine release ##################################################
        if bg_params['plastic_syn']:
            DA_spikes  =nest.Create('poisson_generator')
            my_parrot = nest.Create('parrot_neuron')
            nest.Connect(DA_spikes,my_parrot)
            nest.Connect(my_parrot,bg_params['vt']) #syn_weight=1.

        ############### simulation time parameters ############################
        S1_lag = 300.  #PTN
        M1_lag = 300.
        S1_M1_overlap = 200.
        dop_lag = 10.
        to_rew_lag = 500. 
        from_rew_lag = 500.

        my_spikes_read_shift = 100.#30. #20.

        ###############   map the stimulus as needed and runs M episodes #####################################
        dt_task = S1_lag + M1_lag - S1_M1_overlap #(ms)　
        dt_reward = to_rew_lag + dop_lag + from_rew_lag #dt_task_no_stim = 100.
        start_point = sim_params['initial_ignore'] + dt_task + 100. #(+100. so it starts from 1000.)
        trials_nb = sim_params['sim_model'][sim_model_on]['trials_nb']
        nest.Simulate(start_point) #500ms (initial ignore) the network is silent due not input in CTX by PG.

        ############## Choice reaching ################################################
        print('circle centers: ',bg_params['circle_center']) # print all the channels
        xtarget = np.array(bg_params['circle_center'][1:4])  #select the top 3 targets
        print('targets involve: ',xtarget) 
        target = 1  # target index
        print('prefered target: ',xtarget[1]) #the preferred target, Xd

        ######################################
        ### Initialization (from task) #######
        ######################################
        pm = task.pmm2(mass=0.1,friction=0.01)#03) #mass=0.1, friction=0.2
        ql = task.QLa(len(xtarget))  # 3 actions
        #N = trials_nb  # trials
        a = np.zeros(trials_nb, dtype=int)
        a_GPi = np.zeros(trials_nb, dtype=int)
        r = np.zeros(trials_nb)
        q = np.zeros((trials_nb, ql.na))
        dopa_flag = False
        MSN_fr, GPi_fr, TH_M1_IZ_TC_fr = [],[], []
        push_xy, next_pos, veloc = [],[],[]
        current_pos = [[0.,0.]] + bg_params['circle_center']#np.array([[0.,0.],[0.,0.05],[0.,0.1],[0.,0.15],[0.,0.2],[0.,0.23]]) #current position x
        pm.x = np.array(current_pos[0])
        pm.v=np.array([0.,0.])

        

        #Create files to track activity of each pop:
        with open('./log/'+'MSN_fr.txt', 'a') as file:
            file.write('channel_1(R) channel_2(C) channel_3(L)'+'\n')
        with open('./log/'+'GPi_fr.txt', 'a') as file:
            file.write('channel_1(R) channel_2(C) channel_3(L)'+'\n')
        with open('./log/'+'TH_M1_IZ_TC_fr.txt', 'a') as file:
            file.write('channel_1(R) channel_2(C) channel_3(L)'+'\n')
        with open('./log/'+'push_xy.txt', 'a') as file:
            file.write('x y'+'\n')
        with open('./log/'+'next_pos.txt', 'a') as file:
            file.write('x y'+'\n')
        with open('./log/'+'q_val.txt', 'a') as file:
            file.write('channel_1(R) channel_2(C) channel_3(L)'+'\n')
        with open('./log/'+'reward.txt', 'a') as file:
            file.write('Yes(1.5)-No(0)'+'\n')
        with open('./log/'+'action_sel.txt', 'a') as file:
            file.write('SelectedChannel(1-R,2-C,3-L)'+'\n')
        with open('./log/'+'action_sel_GPi.txt', 'a') as file:
            file.write('SelectedChannel(1-R,2-C,3-L)'+'\n')
        with open('./log/'+'velocity.txt', 'a') as file:
            file.write('x y'+'\n')
        

        ####3
        env = gym.make("MountainCar-v0")
        goal_distance = []
        #####
        for i in np.arange(trials_nb):

            observation = env.reset() 

            #SNN model
            #direction ='CU' 
            print('current_pos: ',current_pos[0])
            
            S1_start = start_point + (S1_lag + M1_lag - S1_M1_overlap + to_rew_lag + dop_lag + from_rew_lag) * float(i) #start_point + i*dt_task +  dt_task_no_stim  #start_point+(new_delta_t*trial_counter)
            S1_end = S1_start + S1_lag #start_point + (i+1)*dt_task - dt_task_no_stim  #start_point+S1_span+(new_delta_t*trial_counter)
            M1_start = S1_start + S1_lag - S1_M1_overlap #start_point + i*dt_task +  dt_task_no_stim   #start_point+S1_span-S1M1_overlap+(new_delta_t*trial_counter)
            M1_end = M1_start + M1_lag #start_point + (i+1)*dt_task - dt_task_no_stim    #start_point+S1_span+M1_span-S1M1_overlap+(new_delta_t*trial_counter)
            dopa_start = M1_end + to_rew_lag #start_point + (S1_lag + M1_lag - S1_M1_overlap + to_rew_lag) * float(i) #S1_start #start_point + i*dt_task                  #new_delta_t-dopa_span/2.+(new_delta_t*trial_counter)
            dopa_end = dopa_start + dop_lag #S1_start + 5.#start_point + (i+1)*dt_task                 #dopa_start + 10.

           
            ###### appply the stimulus to the channels , M2 ####################### 
            print('M2 start: ',S1_start,'  M2 end: ',S1_end,'\n')
            #_ = nest_routine.set_stimulus(current_position=current_pos[0], time_start=S1_start, time_stop=S1_end, psg=S1_L5A_psg, stimulus=stimulus_L5A,centers=bg_params['circle_center']) #S1
            #_ = nest_routine.set_stimulus(current_position=current_pos[0], time_start=S1_start, time_stop=S1_end, psg=S1_L5B_psg, stimulus=stimulus_L5B,centers=bg_params['circle_center']) #S1 
            _ = nest_routine.set_stimulus(current_position=current_pos[0], time_start=S1_start, time_stop=S1_end, psg=CSN_psg, stimulus=stimulus_CSN,centers=bg_params['circle_center']) #M2           

            #print('status: S1_L5B_psg : ',nest.GetStatus(S1_L5B_psg),'\n')
            #print('status: S1_L5A_psg : ',nest.GetStatus(S1_L5A_psg),'\n')
            #print('status: CSN_psg (M2) : ',nest.GetStatus(CSN_psg),'\n')

            #apply stimulus from S1 to the current hand position neurons ####

            #### calculate distance from current place to reference points:
            print('S1 start: ',S1_start,'    S1 end:',S1_end-S1_M1_overlap)
            nearest_idx = nest_routine.get_center_receptive_field(current_pos[0],grid_centers)
            print('nearest_idx: ', nearest_idx) 
            _ = nest_routine.apply_fixed_stimulus(time_start=S1_start, time_stop=S1_end-S1_M1_overlap, psg=PTN_psg[nearest_idx], stimulus=rate_Poisson) #M1

            #_ = nest_routine.apply_direction_stimulus_generic(direction='CU_MSN_M2', time_start=S1_start, time_stop=S1_end, psg=M2_to_MSN_psg, stimulus=stimulus_M2_to_MSN) #M1
            
            
            ######## simulation ##########################
            nest.Simulate(dt_task)
            #nest.Simulate(dt_task)                       #
            ##############################################

            ####### For MSN (get mean firing rate )####################################
            sp_MSN = nest.GetStatus(detectors['MSN'],'events')[0] # will get data for MSN_columns_gids
            idx_dt_MSN = np.where(np.logical_and(sp_MSN['times']>=S1_start, sp_MSN['times']<=S1_end))
            print('MSN')
            MSN_fr.append(nest_routine.get_targets_mean_rates(MSN_columns_gids,sp_MSN,idx_dt_MSN,S1_start,S1_end))
            ####### For GPi (get mean firing rate ) ####################################
            sp_GPi = nest.GetStatus(detectors['GPi'],'events')[0] # will get data for MSN_columns_gids
            idx_dt_GPi = np.where(np.logical_and(sp_GPi['times']>=S1_start+my_spikes_read_shift, sp_GPi['times']<=S1_end+my_spikes_read_shift))  #M1
            print('GPi')
            GPi_fr.append(nest_routine.get_targets_mean_rates(GPi_columns_gids,sp_GPi,idx_dt_GPi,S1_start+my_spikes_read_shift,S1_end+my_spikes_read_shift)) #M1
            ###### Get TH mean firing rates and calculate next push ####################
            sp_TH_M1_IZ_TC = nest.GetStatus(detectors['TH_M1_IZ_thalamic_nucleus_TC'],'events')[0] # will get data for MSN_columns_gids
            idx_dt_TH_M1_IZ_TC = np.where(np.logical_and(sp_TH_M1_IZ_TC['times']>=S1_start+my_spikes_read_shift, sp_TH_M1_IZ_TC['times']<=S1_end+my_spikes_read_shift)) #M1
            print('TH')
            TH_M1_IZ_TC_fr.append(nest_routine.get_targets_mean_rates(TH_M1_IZ_TC_columns_gids,sp_TH_M1_IZ_TC,idx_dt_TH_M1_IZ_TC,S1_start+my_spikes_read_shift,S1_end+my_spikes_read_shift)) #M1

            
            #push_xy.append(nest_routine.get_next_pos(xtarget,TH_M1_IZ_TC_fr_max,current_position=current_pos[0],delta_x=0.1))#delta_x=0.1))
            push_xy.append(nest_routine.get_next_pos(xtarget,TH_M1_IZ_TC_fr[-1],current_position=current_pos[0],delta_x=0.1))#delta_x=0.1))


            print('TH_M1_IZ_TC_fr[-1]: ',TH_M1_IZ_TC_fr[-1])
            print('GPi_fr[-1]: ',GPi_fr[-1])
            print('push_xy[-1]: ',push_xy[-1])
            
            actions_idx = np.argsort(GPi_fr[-1])

            
            if actions_idx[0]==0:
                perform_this = [2] #move right #[2,1] #right-center
            if actions_idx[0]==1:
                perform_this = [1] #dont move at all #[0,1] #left-center
            if actions_idx[0]==2:
                perform_this = [0] #move to left #[0,2] #left-right

            
            ######### open gym ai ###########
            print('observation:  ',observation)
            
            accumulator = 0
            for _ in range(100):
                for kk in perform_this: #actions_idx[:2]: #only 2 actions are ptomoted , 0: left and 2 rights.
                    env.render()                                                          
                    observation, reward, done, info = env.step(kk)  
                    accumulator = accumulator + abs(observation[0]-0.6)
                    print("=" * 10)
                    print("action=",kk)
                    print("observation=",observation)
                    print("reward=",reward)
                    print("done=",done)
                    print("info=",info)
            goal_distance.append(accumulator)
                
            print('goal_distance:    ',goal_distance[-1])
            ###################################################################################
            ql.Q = np.array(MSN_fr[-1]) #last values of MSN firing rates.

            a[i] = ql.select()  #selection by multinomial based on MSN (Q-value)
            a_GPi[i] =  np.argmin(GPi_fr[-1]) #selection by GPi  #also test the selection by TH.

            r[i] = task.reward(a_GPi[i]) #task.reward(a[i])
            #### update ##########
            if actions_idx[0]==2 and i>0:#actions_idx[0]==0 and actions_idx[1]==2: #goal_distance[-1]<=10.: #0.1 #r[i]==1.5:
            #if actions_idx[0]==0:    
                dopa_flag = True
                print('reward ....')
                dopa_rate = 80.#130. #150.
                #if actions_idx[1]==2:
                #    dopa_rate = dopa_rate + 60.
                print('dopa_rate: ',dopa_rate)
            else:
                dopa_flag = False
                print('NO reward .... for episode: '+str(i))
            #######################
            q[i] = ql.Q
            #######################
            
            ######### Chimera end ################################################################

    

            with open('./log/'+'MSN_fr.txt', 'a') as file:
                file.write(str(MSN_fr[-1][0])+' '+str(MSN_fr[-1][1])+' '+str(MSN_fr[-1][2])+'\n')
            with open('./log/'+'GPi_fr.txt', 'a') as file:
                file.write(str(GPi_fr[-1][0])+' '+str(GPi_fr[-1][1])+' '+str(GPi_fr[-1][2])+'\n')
            with open('./log/'+'TH_M1_IZ_TC_fr.txt', 'a') as file:
                file.write(str(TH_M1_IZ_TC_fr[-1][0])+' '+str(TH_M1_IZ_TC_fr[-1][1])+' '+str(TH_M1_IZ_TC_fr[-1][2])+'\n')
            with open('./log/'+'q_val.txt', 'a') as file:
                file.write(str(q[i][0])+' '+str(q[i][1])+' '+str(q[i][2])+'\n')
            with open('./log/'+'reward.txt', 'a') as file:
                file.write(str(r[i])+'\n')
            with open('./log/'+'action_sel.txt', 'a') as file:
                file.write(str(a[i])+'\n')
            with open('./log/'+'action_sel_GPi.txt', 'a') as file:
                file.write(str(a_GPi[i])+'\n')

            ########## dopamine-based reward settings ##########################################
            print('dopa start: ',dopa_start,'dopa end: ',dopa_end,'\n')
            if bg_params['plastic_syn'] and dopa_flag:
                _ = nest.SetStatus(DA_spikes,{'rate':dopa_rate,'start':dopa_start,'stop':dopa_end}) #rate: 150., 20 100., 200., 300.
                print('status DA_spikes : ',nest.GetStatus(DA_spikes),'\n')
            ####################################################################################
            nest.Simulate(dt_reward)



        env.close()

        #### save all related files for plotting #####
        np.savetxt('./log/q_val.txt',q)
        np.savetxt('./log/MSN_fr.txt',MSN_fr)
        np.savetxt('./log/reward.txt',r)
        np.savetxt('./log/action_sel.txt',a)
        np.savetxt('./log/action_sel_GPi.txt',a_GPi)
        np.savetxt('./log/GPi_fr.txt',GPi_fr)
        np.savetxt('./log/TH_M1_IZ_TC_fr.txt',TH_M1_IZ_TC_fr)
        np.savetxt('./log/goal_distance.txt',goal_distance)
        


    
    else:
        print ('wrong model set')

    if sim_regions['BG']:
        for layer_name in bg_layers.keys():
            rate = nest_routine.average_fr(detectors[layer_name], sim_params['simDuration']-sim_params['start_time_sp'],nest_routine.count_layer(bg_layers[layer_name]))
            print('Layer ' + layer_name + " fires at " + str(rate) + " Hz")
    if sim_regions['S1']:
        for layer_name in ctx_layers.keys():
            rate = nest_routine.average_fr(detectors[layer_name], sim_params['simDuration'],nest_routine.count_layer(ctx_layers[layer_name]))
            print('Layer ' + layer_name + " fires at " + str(rate) + " Hz")
    if sim_regions['M1']:
        for layer_name in ctx_M1_layers.keys():
            rate = nest_routine.average_fr(detectors[layer_name], sim_params['simDuration'],nest_routine.count_layer(ctx_M1_layers[layer_name]))
            print('Layer ' + layer_name + " fires at " + str(rate) + " Hz")
    if sim_regions['CB_S1']:
        for layer_name in cb_layers_S1.keys():
            rate = nest_routine.average_fr(detectors[layer_name], sim_params['simDuration'],
                                           nest_routine.count_layer(cb_layers_S1[layer_name]))
            print('Layer ' + layer_name + " fires at " + str(rate) + " Hz")
    if sim_regions['CB_M1']:
        for layer_name in cb_layers_M1.keys():
            rate = nest_routine.average_fr(detectors[layer_name], sim_params['simDuration'],
                                           nest_routine.count_layer(cb_layers_M1[layer_name]))
            print('Layer ' + layer_name + " fires at " + str(rate) + " Hz")

    if sim_regions['TH_S1']:
        for layer_name in th_layers['TH_S1_EZ'].keys():
            rate = nest_routine.average_fr(detectors['TH_S1_EZ' + '_' + layer_name], sim_params['simDuration'],
                                           nest_routine.count_layer(th_layers['TH_S1_EZ'][layer_name]))
            print('Layer ' + 'TH_S1_EZ'+layer_name + " fires at " + str(rate) + " Hz")
    if sim_regions['TH_S1']:
        for layer_name in th_layers['TH_S1_IZ'].keys():
            rate = nest_routine.average_fr(detectors['TH_S1_IZ' + '_' + layer_name], sim_params['simDuration'],
                                           nest_routine.count_layer(th_layers['TH_S1_IZ'][layer_name]))
            print('Layer ' + 'TH_S1_IZ'+ layer_name + " fires at " + str(rate) + " Hz")
    if sim_regions['TH_M1']:
        for layer_name in th_layers['TH_M1_EZ'].keys():
            rate = nest_routine.average_fr(detectors['TH_M1_EZ' + '_' + layer_name], sim_params['simDuration'],
                                           nest_routine.count_layer(th_layers['TH_M1_EZ'][layer_name]))
            print('Layer ' + 'TH_M1_EZ' +layer_name + " fires at " + str(rate) + " Hz")
    if sim_regions['TH_M1']:
        for layer_name in th_layers['TH_M1_IZ'].keys():
            rate = nest_routine.average_fr(detectors['TH_M1_IZ' + '_' + layer_name], sim_params['simDuration'],
                                           nest_routine.count_layer(th_layers['TH_M1_IZ'][layer_name]))
            print('Layer ' + 'TH_M1_IZ' + layer_name + " fires at " + str(rate) + " Hz")

if __name__ == '__main__':
    main()



