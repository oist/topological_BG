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
import json
import random



def main():
    # 1) reads parameters
    sim_params = fetch_params.read_sim()
    bg_params = fetch_params.read_bg()

    # 1.5) initialize nest
    nest_routine.initialize_nest(sim_params)

    # 2) instantiates regions
    if sim_regions['BG']:
        if sim_params['channels']: #if True
            bg_params['channels'] = True #for channels input tasks
        else:
            bg_params['channels'] = False #resting state
        bg_params['circle_center'] = nest_routine.get_channel_centers(sim_params, hex_center=[0, 0],
                                                                ci=sim_params['channels_nb'],
                                                                hex_radius=sim_params['hex_radius'])
        bg_layers,ctx_bg_input = ini_all.instantiate_bg(bg_params, fake_inputs=True,
                                                                ctx_inputs=None,scalefactor=sim_params['scalefactor'])

    

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



    if sim_model_on=='resting_state':

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
        #select a subset of 200 neurons for each channel
        CSN_columns_gids_subset = []
        for k1 in CSN_columns_gids:
            CSN_columns_gids_subset.append(random.sample(k1,sample_input_size))
        _ = nest_routine.connect_psg_to_channels(CSN_columns_gids_subset, CSN_psg, CSN_syn)

        PTN_columns_gids = nest_routine.get_columns_data('PTN', bg_params['circle_center'],sim_params['channels_radius']) 
        PTN_psg, PTN_syn = nest_routine.create_psg_channels(syn_weight=1., syn_delay=2.,channels_nb=sim_params['channels_nb']) #0.3 #tested: 27.5 (looks caothic), 7.5, 3., 1., 0.5, 0.1, 5.#original: 11.5
        #select a subset of 200 neurons for each channel
        PTN_columns_gids_subset = []
        for k2 in PTN_columns_gids:
            PTN_columns_gids_subset.append(random.sample(k2,sample_input_size))
        _ = nest_routine.connect_psg_to_channels(PTN_columns_gids_subset, PTN_psg, PTN_syn)

        MSN_columns_gids = nest_routine.get_columns_data('MSN', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN
        MSN_d1_columns_gids = nest_routine.get_columns_data('MSN_d1', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN_d1
        MSN_d2_columns_gids = nest_routine.get_columns_data('MSN_d2', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN_d2

        GPi_columns_gids = nest_routine.get_columns_data('GPi', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of GPi

        channel_1_idx = 1
        channel_2_idx = 2 

        nest.Simulate(1000.)
        stim_dt = 500.
        stim_break = 500.
        episode = 0
        MSN_fr, MSN_d1_fr, MSN_d2_fr, GPi_fr, PTN_fr, CSN_fr = [],[],[],[], [], []
        PTN_CSN_ratio = 1.8

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

        for i in np.arange(11):  #channel 1
            for j in np.arange(11):   #channel 2

                Poisson_ch1 = i*1.7 
                Poisson_ch2 = j*1.7 
                
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

                ####### For MSN (get the mean firing rate )####################################
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
        print('finished')
   
    elif sim_model_on == 'plasticity':
        if bg_params['plastic_syn']:

            DA_spikes_d1  =nest.Create('poisson_generator')
            my_parrot_d1 = nest.Create('parrot_neuron')
            nest.Connect(DA_spikes_d1,my_parrot_d1)
            nest.Connect(my_parrot_d1,bg_params['vt_d1']) #syn_weight=1.

            DA_spikes_d2  =nest.Create('poisson_generator')
            my_parrot_d2 = nest.Create('parrot_neuron')
            nest.Connect(DA_spikes_d2,my_parrot_d2)
            nest.Connect(my_parrot_d2,bg_params['vt_d2']) #syn_weight=1.

            ##### get columns data , create PG and connect ############
            sample_input_size = 300 

            CSN_columns_gids = nest_routine.get_columns_data('CSN', bg_params['circle_center'],sim_params['channels_radius']) 
            CSN_psg, CSN_syn = nest_routine.create_psg_channels(syn_weight=1., syn_delay=2.,channels_nb=sim_params['channels_nb']) #0.3 #tested: 27.5 (looks caothic), 7.5, 3., 1., 0.5, 0.1, 5.#original: 11.5
            #select a subset of 300 neurons for each channel
            CSN_columns_gids_subset = []
            for k1 in CSN_columns_gids:
                CSN_columns_gids_subset.append(random.sample(k1,sample_input_size))
            _ = nest_routine.connect_psg_to_channels(CSN_columns_gids_subset, CSN_psg, CSN_syn)
            
            PTN_columns_gids = nest_routine.get_columns_data('PTN', bg_params['circle_center'],sim_params['channels_radius']) 
            PTN_psg, PTN_syn = nest_routine.create_psg_channels(syn_weight=1., syn_delay=2.,channels_nb=sim_params['channels_nb']) #0.3 #tested: 27.5 (looks caothic), 7.5, 3., 1., 0.5, 0.1, 5.#original: 11.5
            #select a subset of 300 neurons for each channel
            PTN_columns_gids_subset = []
            for k2 in PTN_columns_gids:
                PTN_columns_gids_subset.append(random.sample(k2,sample_input_size))
            _ = nest_routine.connect_psg_to_channels(PTN_columns_gids_subset, PTN_psg, PTN_syn)

            

            
            MSN_columns_gids = nest_routine.get_columns_data('MSN', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN
            MSN_d1_columns_gids = nest_routine.get_columns_data('MSN_d1', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN
            MSN_d2_columns_gids = nest_routine.get_columns_data('MSN_d2', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of MSN

            GPi_columns_gids = nest_routine.get_columns_data('GPi', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of GPi
            
            FSI_columns_gids = nest_routine.get_columns_data('FSI', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of FSI
            GPe_columns_gids = nest_routine.get_columns_data('GPe', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of GPe
            STN_columns_gids = nest_routine.get_columns_data('STN', bg_params['circle_center'],sim_params['channels_radius'])  #this gids are used for calculate mean rate of STN

            channel_1_idx = 1
            channel_2_idx = 2 

            ###### added to track the stimuli #####
            PTN_sd_1 = nest.Create("spike_detector", params={'label': 'PTN_300_ch_1',"withgid": True, "withtime": True, "to_file": True, "fbuffer_size":8192})
            nest.Connect(pre=[item[0] for item in PTN_columns_gids_subset[channel_1_idx]], post=PTN_sd_1)

            PTN_sd_2 = nest.Create("spike_detector", params={'label': 'PTN_300_ch_2',"withgid": True, "withtime": True, "to_file": True, "fbuffer_size":8192})
            nest.Connect(pre=[item[0] for item in PTN_columns_gids_subset[channel_2_idx]], post=PTN_sd_2)

            CSN_sd_1 = nest.Create("spike_detector", params={'label': 'CSN_300_ch_1',"withgid": True, "withtime": True, "to_file": True, "fbuffer_size":8192})
            nest.Connect(pre=[item[0] for item in CSN_columns_gids_subset[channel_1_idx]], post=CSN_sd_1)

            CSN_sd_2 = nest.Create("spike_detector", params={'label': 'CSN_300_ch_2',"withgid": True, "withtime": True, "to_file": True, "fbuffer_size":8192})
            nest.Connect(pre=[item[0] for item in CSN_columns_gids_subset[channel_2_idx]], post=CSN_sd_2)
            ###### ####################################

            global modulated
            modulated = False


            nest.Simulate(1000.)
            stim_dt = 1000.
            stim_break = 500. #time after stimulation offset
            dopa_break = 1000. #time after dopamine release 
            dopa_dt = 300. 
            MSN_fr, GPi_fr, PTN_fr, CSN_fr, MSN_d1_fr, MSN_d2_fr, CMPf_fr, FSI_fr, GPe_fr, STN_fr = [],[],[],[],[],[],[],[],[],[]
            CS_plus, CS_minus = [],[]
           
            Istrong = 8.0 
            Iweak = 6.0 
            Poisson_ch1 = Istrong * 1.7 
            Poisson_ch2 = Iweak * 1.7 
            dopa_rate = 20.0 

            PTN_CSN_ratio = 1.8 
            hits = 0

            current_task = [] #1:CS+, -1:CS-
            con_ep = 100 
            dis_ep= 100 
            gen_ep= 30 
            cond_from = con_ep + dis_ep + gen_ep
            total_episodes = cond_from
            dc = 0.00056 

            for episode in np.arange(total_episodes):
                if episode > cond_from: 
                    if random.random() > 1.0: # this is CS+
                        CS_plus.append(1)
                        CS_minus.append(0)
                        Poisson_ch1 = Istrong * 1.7 
                        Poisson_ch2 = Iweak * 1.7 
                        current_task.append(1)
                    else:  ##this is CS-
                        CS_minus.append(1)
                        CS_plus.append(0)
                        Poisson_ch1 = Iweak * 1.7 
                        Poisson_ch2 = Istrong * 1.7 
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
  
                if episode > ( con_ep + gen_ep):
                    with open('./log/' + 'do_tests.txt', 'r') as file:
                        for lines in file:
                            do_tests = lines
                    print('do_tests: ',do_tests,'length: ',len(do_tests))
                    if do_tests[:-1]=='T2':
                        Poisson_ch1 = Iweak * 1.7 
                        Poisson_ch2 = Istrong * 1.7 
                    if do_tests[:-1]=='T1': 
                        Poisson_ch1 = Istrong * 1.7 
                        Poisson_ch2 = Iweak * 1.7 
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
                        Poisson_ch1 = Istrong * 1.7 
                        Poisson_ch2 = Iweak * 1.7 

                    if do_tests[:-1]=='T' and (episode==(con_ep + 10)  or episode==(con_ep + 20)):
                        do_tests= 'T2 '
                        Poisson_ch1 = Iweak * 1.7 
                        Poisson_ch2 = Istrong * 1.7 
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

                if episode < total_episodes-2 and do_tests[:-1]=='F':  

                    if current_task[-1] == 1: #CS+ task
                        if GPi_fr[-1][channel_1_idx] < GPi_fr[-1][channel_2_idx]: #if channel 1 is winner
                            _ = nest.SetStatus(DA_spikes_d1,{'rate':dopa_rate,'start':dopa_start,'stop':dopa_end}) 
                            print('episode: ',episode,' ,task: ',current_task[-1],', Dopamine burst from vt_d1:',dopa_start,'  to:',dopa_end,'\n')
                
                    if current_task[-1] == -1 and MSN_d2_fr[-1][2]<70.0:  #CS- task 
                        if GPi_fr[-1][channel_1_idx] < GPi_fr[-1][channel_2_idx]: #if channel 1 is winner
                            _ = nest.SetStatus(DA_spikes_d2,{'rate':dopa_rate,'start':dopa_start,'stop':dopa_end+100.}) #+ 100 ms for the case of dip
                            print('episode: ',episode,' ,task: ',current_task[-1],', Dopamine Release from vt_d2:',dopa_start,'  to:',dopa_end+100.,'\n')
                            
                            #modulation of D2 neurons synapses.
                            if not(modulated): #then modulate.
                                my_gain = bg_params['MSND2_mod'] 

                                my_target = nest.GetNodes(bg_layers['MSN'])[0][:] #'MSNs'
                                my_source = nest.GetNodes(bg_layers['MSN'])[0][int(bg_params['nbMSN']/2):] #'MSN_d2'
                                conn = nest.GetConnections(source=my_source, target=my_target)
                                old_weight = nest.GetStatus(conn)[0]['weight']

                                newweight = old_weight*my_gain 
                                nest.SetStatus(conn, {'weight': newweight})

                                my_target_1 = nest.GetNodes(bg_layers['GPe'])[0]
                                conn_1 = nest.GetConnections(source=my_source, target=my_target_1)
                                old_weight_1 = nest.GetStatus(conn_1)[0]['weight']

                                newweight_1 = old_weight_1*my_gain 
                                nest.SetStatus(conn_1, {'weight': newweight_1})

                                my_target_2 = nest.GetNodes(bg_layers['GPi'])[0]
                                conn_2 = nest.GetConnections(source=my_source, target=my_target_2)
                                old_weight_2 = nest.GetStatus(conn_2)[0]['weight']

                                newweight_2 = old_weight_2*my_gain #2.0 #4.0 #2.0
                                nest.SetStatus(conn_2, {'weight': newweight_2})

                                modulated = True
                        else:
                            print('episode: ',episode,' ,task: ',current_task[-1],'no dopamine release, selected channel 2')
                else:
                    print('testing GPi response without Dop burst or dip')


                if episode == con_ep: 
                    cond_from = episode


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

            print('finished')

        else:
            print('error: plastic_syn flag is deactivate at bgParams')

    
    else:
        print ('wrong model set')

    if sim_regions['BG']:
        for layer_name in bg_layers.keys():
            rate = nest_routine.average_fr(detectors[layer_name], sim_params['simDuration']-sim_params['start_time_sp'],nest_routine.count_layer(bg_layers[layer_name]))
            print('Layer ' + layer_name + " fires at " + str(rate) + " Hz")
    
if __name__ == '__main__':
    main()



