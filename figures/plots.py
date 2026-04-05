import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import time

# Path to the reference simulation log with connectivity files
DATA_PATH = './2021_05_03_16_58_22_xp000000/log_focused_MSN/'


def hex_corner(center, size, i):
    angle_deg = 60 * i - 30
    angle_rad = np.pi / 180 * angle_deg
    return [center[0] + size * np.cos(angle_rad), center[1] + size * np.sin(angle_rad)]


def draw_hex_grid(ax, radius_circle=0.12):
    """Draw the hexagonal channel grid and channel circles."""
    for k in [[0,0],[0.240*2.-0.06433,0],[0.-(0.240*2.-0.06433),0],
              [0.240*1.-0.06433/2.,0.240*2-0.12],[0.240*1.-0.06433/2.,-(0.240*2-0.12)],
              [-(0.240*1.-0.06433/2.),0.240*2-0.12],[-(0.240*1.-0.06433/2.),-(0.240*2-0.12)]]:
        my_centers = []
        for i in np.arange(6):
            x_y = hex_corner(k, 0.240, i)
            ax.plot(x_y[0], x_y[1], '.', color='k')
            my_centers.append(x_y)
            if k == [0, 0]:
                my_circle = plt.Circle(x_y, radius_circle, color='k', alpha=0.5,
                                       linewidth=2, fill=False)
                ax.add_patch(my_circle)
        ax.plot(np.array(my_centers)[:,0], np.array(my_centers)[:,1], '-k', linewidth=0.4)
        ax.plot([my_centers[-1][0], my_centers[0][0]],
                [my_centers[-1][1], my_centers[0][1]], '-k', linewidth=0.4)


def Figure_1_E(save=True):
    """
    Figure 1E: Cortical afferents from channels 1 and 2 (CSN + PTN)
    projecting onto MSN_d1 and MSN_d2 neurons in the striatum.
    """
    fig, ax = plt.subplots(figsize=(8., 8.))

    MSN_d1_pos = np.loadtxt(DATA_PATH + 'MSN_d1.txt')
    MSN_d2_pos = np.loadtxt(DATA_PATH + 'MSN_d2.txt')
    CSN_pos    = np.loadtxt(DATA_PATH + 'CSN.txt')
    PTN_pos    = np.loadtxt(DATA_PATH + 'PTN.txt')

    for i, c, a in zip([1, 2], ['b', 'crimson'], [0.2, 0.2]):
        CSN_d1 = np.loadtxt(DATA_PATH + 'CSN_c'+str(i)+'_to_MSN_d1.txt')
        CSN_d2 = np.loadtxt(DATA_PATH + 'CSN_c'+str(i)+'_to_MSN_d2.txt')
        PTN_d1 = np.loadtxt(DATA_PATH + 'PTN_c'+str(i)+'_to_MSN_d1.txt')
        PTN_d2 = np.loadtxt(DATA_PATH + 'PTN_c'+str(i)+'_to_MSN_d2.txt')

        my_idx_2 = np.isin(MSN_d1_pos[:,0], CSN_d1[:,0])
        my_idx_4 = np.isin(MSN_d2_pos[:,0], CSN_d2[:,0])
        my_idx_6 = np.isin(MSN_d1_pos[:,0], PTN_d1[:,0])
        my_idx_8 = np.isin(MSN_d2_pos[:,0], PTN_d2[:,0])

        ax.scatter(MSN_d1_pos[:,1][my_idx_2], MSN_d1_pos[:,2][my_idx_2], marker='o', s=28, color=c, alpha=a)
        ax.scatter(MSN_d2_pos[:,1][my_idx_4], MSN_d2_pos[:,2][my_idx_4], marker='x', s=12, color=c, alpha=0.4)
        ax.scatter(MSN_d1_pos[:,1][my_idx_6], MSN_d1_pos[:,2][my_idx_6], marker='o', s=28, color=c, alpha=a)
        ax.scatter(MSN_d2_pos[:,1][my_idx_8], MSN_d2_pos[:,2][my_idx_8], marker='x', s=12, color=c, alpha=0.4)

    draw_hex_grid(ax)
    ax.set_xlim([-0.5, 0.5])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticks([])
    ax.set_yticks([])

    legend_elements = [
        Line2D([0],[0], ls='', marker='o', color='crimson', markerfacecolor='crimson', markersize=5, alpha=0.4, label='Ch2 cortical afferents'),
        Line2D([0],[0], ls='', marker='o', color='b',       markerfacecolor='b',       markersize=5, alpha=0.4, label='Ch1 cortical afferents'),
        Line2D([0],[0], ls='', marker='o', color='k',       markerfacecolor='white',   markersize=5, alpha=0.5, label='D1-MSN'),
        Line2D([0],[0], ls='', marker='x', color='k',                                  markersize=5, alpha=0.5, label='D2-MSN'),
    ]
    fig.legend(handles=legend_elements, loc=(0.22, 0.88), labelspacing=0.1, ncol=2, fontsize=14, frameon=False)

    if save:
        plt.savefig('fig_E_'+str(time.time())+'.png', dpi=300)
    plt.show()


def Figure_1_F(save=True):
    """
    Figure 1F: Structural asymmetry from a single MSN_D1 presynaptic neuron (square).
    Shows postsynaptic MSN_d1 targets (circles) and MSN_d2 targets (crosses)
    for channels 1 and 2 — illustrating that MSN_D1->MSN_D2 connections are sparse.
    """
    fig, ax = plt.subplots(figsize=(8., 8.))

    MSN_d1_pos = np.loadtxt(DATA_PATH + 'MSN_d1.txt')
    MSN_d2_pos = np.loadtxt(DATA_PATH + 'MSN_d2.txt')

    for i, c, d1 in zip([1, 2], ['b', 'crimson'], [150, 105]):
        d1_d1 = np.loadtxt(DATA_PATH + 'MSN_d1_c'+str(i)+'_to_MSN_d1.txt')
        d1_d2 = np.loadtxt(DATA_PATH + 'MSN_d1_c'+str(i)+'_to_MSN_d2.txt')

        # Pick one MSN_d1 presynaptic source neuron
        d1_source = d1_d1[d1, 1]
        idx_source = np.where(MSN_d1_pos[:,0] == d1_source)[0]
        ax.scatter(MSN_d1_pos[:,1][idx_source], MSN_d1_pos[:,2][idx_source],
                   marker='s', s=108, color=c, edgecolor='k', alpha=0.9)

        # Its MSN_d1 postsynaptic targets
        idx_targets = np.where(d1_d1[:,1] == d1_source)[0]
        my_idx_0 = np.isin(MSN_d1_pos[:,0], d1_d1[:,0][idx_targets])
        ax.scatter(MSN_d1_pos[:,1][my_idx_0], MSN_d1_pos[:,2][my_idx_0],
                   marker='o', s=88, color=c, alpha=0.4)

        # Its MSN_d2 postsynaptic targets
        idx_targets = np.where(d1_d2[:,1] == d1_source)[0]
        my_idx_0 = np.isin(MSN_d2_pos[:,0], d1_d2[:,0][idx_targets])
        ax.scatter(MSN_d2_pos[:,1][my_idx_0], MSN_d2_pos[:,2][my_idx_0],
                   marker='x', s=70, color=c, alpha=0.5)

    draw_hex_grid(ax)
    ax.set_xlim([-0.25, 0.45])
    ax.set_ylim([-0.2, 0.5])
    ax.set_xticks([])
    ax.set_yticks([])

    legend_elements = [
        Line2D([0],[0], ls='', marker='o', color='crimson', markerfacecolor='crimson', markersize=5, alpha=0.4, label='Channel 2'),
        Line2D([0],[0], ls='', marker='o', color='b',       markerfacecolor='b',       markersize=5, alpha=0.4, label='Channel 1'),
        Line2D([0],[0], ls='', marker='o', color='k',       markerfacecolor='white',   markersize=5, alpha=0.5, label="Postsynaptic D1-MSN's"),
        Line2D([0],[0], ls='', marker='x', color='k',                                  markersize=5, alpha=0.5, label="Postsynaptic D2-MSN's"),
        Line2D([0],[0], ls='', marker='s', color='k',       markerfacecolor='white',   markersize=8, alpha=0.5, label='Presynaptic D1-MSN'),
    ]
    fig.legend(handles=legend_elements, loc=(0.12, 0.88), labelspacing=0.1, ncol=2, fontsize=14, frameon=False)

    if save:
        plt.savefig('fig_F_'+str(time.time())+'.png', dpi=300)
    plt.show()


def Figure_1_G(save=True):
    """
    Figure 1G: Structural asymmetry from a single MSN_D2 presynaptic neuron (triangle).
    Shows postsynaptic MSN_d1 targets (circles) and MSN_d2 targets (crosses)
    for channels 1 and 2 — illustrating that MSN_D2 connections are more symmetric.
    """
    fig, ax = plt.subplots(figsize=(8., 8.))

    MSN_d1_pos = np.loadtxt(DATA_PATH + 'MSN_d1.txt')
    MSN_d2_pos = np.loadtxt(DATA_PATH + 'MSN_d2.txt')

    for i, c, d2 in zip([1, 2], ['b', 'crimson'], [150, 105]):
        d2_d1 = np.loadtxt(DATA_PATH + 'MSN_d2_c'+str(i)+'_to_MSN_d1.txt')
        d2_d2 = np.loadtxt(DATA_PATH + 'MSN_d2_c'+str(i)+'_to_MSN_d2.txt')

        # Pick one MSN_d2 presynaptic source neuron
        d2_source = d2_d1[d2, 1]
        idx_source = np.where(MSN_d2_pos[:,0] == d2_source)[0]
        ax.scatter(MSN_d2_pos[:,1][idx_source], MSN_d2_pos[:,2][idx_source],
                   marker='^', s=108, color=c, edgecolor='k', alpha=0.9)

        # Its MSN_d1 postsynaptic targets
        idx_targets = np.where(d2_d1[:,1] == d2_source)[0]
        my_idx_0 = np.isin(MSN_d1_pos[:,0], d2_d1[:,0][idx_targets])
        ax.scatter(MSN_d1_pos[:,1][my_idx_0], MSN_d1_pos[:,2][my_idx_0],
                   marker='o', s=88, color=c, alpha=0.4)

        # Its MSN_d2 postsynaptic targets
        idx_targets = np.where(d2_d2[:,1] == d2_source)[0]
        my_idx_0 = np.isin(MSN_d2_pos[:,0], d2_d2[:,0][idx_targets])
        ax.scatter(MSN_d2_pos[:,1][my_idx_0], MSN_d2_pos[:,2][my_idx_0],
                   marker='x', s=70, color=c, alpha=0.5)

    draw_hex_grid(ax)
    ax.set_xlim([-0.25, 0.45])
    ax.set_ylim([-0.2, 0.5])
    ax.set_xticks([])
    ax.set_yticks([])

    legend_elements = [
        Line2D([0],[0], ls='', marker='o', color='crimson', markerfacecolor='crimson', markersize=5, alpha=0.4, label='Channel 2'),
        Line2D([0],[0], ls='', marker='o', color='b',       markerfacecolor='b',       markersize=5, alpha=0.4, label='Channel 1'),
        Line2D([0],[0], ls='', marker='o', color='k',       markerfacecolor='white',   markersize=5, alpha=0.5, label="Postsynaptic D1-MSN's"),
        Line2D([0],[0], ls='', marker='x', color='k',                                  markersize=5, alpha=0.5, label="Postsynaptic D2-MSN's"),
        Line2D([0],[0], ls='', marker='^', color='k',       markerfacecolor='white',   markersize=8, alpha=0.5, label='Presynaptic D2-MSN'),
    ]
    fig.legend(handles=legend_elements, loc=(0.12, 0.88), labelspacing=0.1, ncol=2, fontsize=14, frameon=False)

    if save:
        plt.savefig('fig_G_'+str(time.time())+'.png', dpi=300)
    plt.show()


def Figure_3_C(save=True):
    """
    Figure C: Generalization and discrimination learning trajectory.
    Top panel: D1-MSN and D2-MSN mean firing rates for Ch1 and Ch2 across episodes.
    Bottom panel: GPi mean firing rates for Ch1 and Ch2 across episodes.

    Three phases separated by vertical dashed lines:
      0  .. 100 : CS+ conditioning  (Ch1 stimulated, DA bursts)
      101 .. 130 : generalization probing (no stimulus; T1/T2 probes at eps +5,+10,+15,+20)
      131 .. 229 : CS- discrimination (Ch2 stimulated, DA dips)

    Data: run4_k_4_l_02_m_1/log1  (kappa=4.0, lambda=0.2, modulation=1.0, no eta modulation)
    """
    DATA_C = './run4_k_4_l_02_m_1/log1/'
    # Phase boundaries
    gen_interv = [100, 30, 100]   # [CS+ eps, generalization eps, CS- eps]
    epis       = 228              # episodes plotted (0..227)
    ch1, ch2   = 1, 2
    msize      = 8

    MSN_d1_fr = np.loadtxt(DATA_C + 'MSN_d1_fr.txt')
    MSN_d2_fr = np.loadtxt(DATA_C + 'MSN_d2_fr.txt')
    GPi_fr    = np.loadtxt(DATA_C + 'GPi_fr.txt')

    params    = np.loadtxt(DATA_C + 'plot_labels.txt', dtype=str)
    my_lambda = params[0].split(',')[1]
    my_kappa  = params[1].split(',')[1]

    eps = np.arange(epis)

    fig, ax = plt.subplots(2, figsize=(12, 12))

    # --- MSN panel ---
    color_b = 'steelblue'
    color_r = 'crimson'

    ax[0].plot(eps, MSN_d1_fr[:epis, ch1], 'o', label='D1-MSN (Ch1), $\\lambda$='+my_lambda,
               alpha=0.7, markeredgewidth=1.1, markeredgecolor=color_b, markersize=msize, markerfacecolor='None')
    ax[0].plot(eps, MSN_d2_fr[:epis, ch1], 'o', label='D2-MSN (Ch1), $\\lambda$='+my_lambda,
               alpha=0.7, markeredgewidth=1.1, markeredgecolor=color_r, markersize=msize, markerfacecolor='None')
    ax[0].plot(eps, MSN_d2_fr[:epis, ch2], '+', label='D2-MSN (Ch2), $\\lambda$='+my_lambda,
               alpha=0.7, markeredgewidth=1.1, markeredgecolor=color_r, markersize=msize, markerfacecolor='None')
    ax[0].plot(eps, MSN_d1_fr[:epis, ch2], '+', label='D1-MSN (Ch2), $\\lambda$='+my_lambda,
               alpha=0.7, markeredgewidth=1.1, markeredgecolor=color_b, markersize=msize, markerfacecolor='None')

    # Generalization test markers (stim A: eps +5,+15; stim B: eps +10,+20)
    g0 = gen_interv[0]
    ax[0].plot([g0+5,  g0+15], [MSN_d1_fr[g0+5,  ch1], MSN_d1_fr[g0+15, ch1]], '*', color='k', label='Gen. Test (stim. A)', ms=msize, alpha=0.5)
    ax[0].plot([g0+10, g0+20], [MSN_d1_fr[g0+10, ch1], MSN_d1_fr[g0+20, ch1]], 's', color='k', label='Gen. Test (stim. B)', ms=msize, alpha=0.5)
    ax[0].plot([g0+5,  g0+15], [MSN_d2_fr[g0+5,  ch2], MSN_d2_fr[g0+15, ch2]], '*', color='k', ms=msize, alpha=0.5)
    ax[0].plot([g0+10, g0+20], [MSN_d2_fr[g0+10, ch2], MSN_d2_fr[g0+20, ch2]], 's', color='k', ms=msize, alpha=0.5)

    # --- GPi panel ---
    ax[1].plot(eps, GPi_fr[:epis, ch1], 'o', label='GPi (Ch1)',
               alpha=0.7, markeredgewidth=1.1, markeredgecolor=color_b, markersize=msize, markerfacecolor='None')
    ax[1].plot(eps, GPi_fr[:epis, ch2], '+', label='GPi (Ch2)',
               alpha=0.7, markeredgewidth=1.1, markeredgecolor=color_r, markersize=msize, markerfacecolor='None')

    ax[1].plot([g0+5,  g0+15], [GPi_fr[g0+5,  ch2], GPi_fr[g0+15, ch2]], '*', color='k', label='Gen. Test (stim. A)', ms=msize, alpha=0.5)
    ax[1].plot([g0+10, g0+20], [GPi_fr[g0+10, ch2], GPi_fr[g0+20, ch2]], 's', color='k', label='Gen. Test (stim. B)', ms=msize, alpha=0.5)
    ax[1].plot([g0+5,  g0+15], [GPi_fr[g0+5,  ch1], GPi_fr[g0+15, ch1]], '*', color='k', ms=msize, alpha=0.5)
    ax[1].plot([g0+10, g0+20], [GPi_fr[g0+10, ch1], GPi_fr[g0+20, ch1]], 's', color='k', ms=msize, alpha=0.5)

    # Phase boundary lines
    for bound in [g0, g0+gen_interv[1], g0+gen_interv[1]+gen_interv[2]]:
        ax[0].axvline(bound, color='k', linestyle='-.', alpha=0.5)
        ax[1].axvline(bound, color='k', linestyle='-.', alpha=0.5)

    ax[0].set_ylim([-0.2, 7.55])
    ax[1].set_ylim([43, 92])

    ax[0].set_title('MSN\n\n$\\kappa$='+my_kappa, fontsize=18)
    ax[1].set_title('GPi', fontsize=18)
    ax[0].set_ylabel('Mean firing rate [Hz]', fontsize=16)
    ax[1].set_ylabel('Mean firing rate [Hz]', fontsize=16)
    ax[1].set_xlabel('Episodes', fontsize=16)

    ax[0].legend(loc='upper left', fontsize=16)
    ax[1].legend(loc='lower left', fontsize=16)

    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)

    ax[0].tick_params(axis='both', which='major', labelsize=14)
    ax[1].tick_params(axis='both', which='major', labelsize=14)

    plt.subplots_adjust(hspace=0.35)
    if save:
        plt.savefig('fig_C_'+str(time.time())+'.png', format='png', dpi=300)
    plt.show()


def Figure_3_D(save=True):
    """
    Figure D: STDP functions for D1-MSN and D2-MSN synapses.
    2x2 grid: (D1 burst, D1 dip) top row; (D2 burst, D2 dip) bottom row.

    STDP curve:
      Δt > 0 (pre before post): A_plus  * exp(-Δt / tau_plus)
      Δt < 0 (post before pre): -A_minus * exp( Δt / tau_minus)

    Parameters (dc = 0.00056, tau = 20 ms):
      D1 DA burst: A_plus=dc,      A_minus=dc/4
      D1 DA dip:   A_plus=-dc/40,  A_minus=dc/40
      D2 DA burst: A_plus=dc/3,    A_minus=dc*3/4
      D2 DA dip:   A_plus=dc,      A_minus=-dc
    """
    dc  = 0.00056
    tau = 20.0
    dt  = np.arange(0., 40., 0.01)

    conditions = [
        # (title,               A_plus,     A_minus,     color,       row, col)
        ('D1-MSN (DA burst)',   dc,          dc/4.0,     'steelblue', 0,   0),
        ('D1-MSN (DA dip)',    -dc/40.0,    dc/40.0,     'steelblue', 0,   1),
        ('D2-MSN (DA burst)',   dc/3.0,      dc*3.0/4.0, 'crimson',   1,   0),
        ('D2-MSN (DA dip)',     dc,         -dc,          'crimson',   1,   1),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(8, 8), sharex=True)

    for title, Aplus, Aminus, color, row, col in conditions:
        ax = axes[row][col]

        # Pre before post (Δt > 0)
        ax.plot(dt,      Aplus  * np.exp(-dt / tau), color=color, linewidth=3)
        # Post before pre (Δt < 0)
        ax.plot(-1.*dt, -Aminus * np.exp(-dt / tau), color=color, linewidth=3)

        # Dotted crosshairs
        ax.plot([0., 0.], [-dc, dc], 'k:')
        ax.plot([-40, 40], [0., 0.], 'k:')

        ax.set_title(title, fontsize=18)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks([0])
        ax.tick_params(axis='both', which='major', labelsize=15)

        # A+ / A- labels on left column only
        if col == 0:
            ax.text(-42, dc * 0.6,  'A+', fontsize=16, va='center')
            ax.text(-42, -dc * 0.6, 'A-', fontsize=16, va='center')

        if row == 1:
            ax.set_xlabel('$\\Delta t$ (ms)', fontsize=17)

    plt.tight_layout()
    if save:
        plt.savefig('fig_D_'+str(time.time())+'.png', format='png', dpi=300)
    plt.show()


def Figure_4_A(save=True):
    """
    Figure A: CS+ conditioning — individual runs + mean±std for D1-MSN (top)
    and GPi (bottom), episodes 0–100.
    Data: ./run4_k_4_l_002_m_1/ (κ=4, λ≈0.0, η=1)
    """
    base_path = './run4_k_4_l_002_m_1/'
    my_color = 'midnightblue'
    from_here = 0
    to_here = 100
    c1 = 1
    c2 = 2

    fig, ax = plt.subplots(2, figsize=(8, 10))

    app_MSN_d1_c1 = np.zeros((to_here - from_here, 1))
    app_MSN_d1_c2 = np.zeros((to_here - from_here, 1))
    app_GPi_c1    = np.zeros((to_here - from_here, 1))
    app_GPi_c2    = np.zeros((to_here - from_here, 1))

    my_lambda = my_kappa = ''
    my_count = 1

    for m in np.arange(20):
        my_path = base_path + 'log' + str(m) + '/'
        try:
            MSN_d1_fr = np.loadtxt(my_path + 'MSN_d1_fr.txt')
            GPi_fr    = np.loadtxt(my_path + 'GPi_fr.txt')

            app_MSN_d1_c1 = np.append(app_MSN_d1_c1,
                MSN_d1_fr[from_here:to_here, c1].reshape(to_here - from_here, 1), axis=1)
            app_MSN_d1_c2 = np.append(app_MSN_d1_c2,
                MSN_d1_fr[from_here:to_here, c2].reshape(to_here - from_here, 1), axis=1)
            app_GPi_c1 = np.append(app_GPi_c1,
                GPi_fr[from_here:to_here, c1].reshape(to_here - from_here, 1), axis=1)
            app_GPi_c2 = np.append(app_GPi_c2,
                GPi_fr[from_here:to_here, c2].reshape(to_here - from_here, 1), axis=1)

            ax[0].plot(np.arange(from_here, to_here), MSN_d1_fr[from_here:to_here, c1],
                       'o', label='Ch1, run-' + str(my_count), alpha=0.6)
            ax[1].plot(np.arange(from_here, to_here), GPi_fr[from_here:to_here, c1],
                       'o', alpha=0.6)

            params = np.loadtxt(my_path + 'plot_labels.txt', dtype=str)
            my_lambda = params[0].split(',')[1]
            my_kappa  = params[1].split(',')[1]
            my_count += 1
        except Exception:
            pass

    # Remove initialisation column
    app_MSN_d1_c1 = np.delete(app_MSN_d1_c1, 0, 1)
    app_MSN_d1_c2 = np.delete(app_MSN_d1_c2, 0, 1)
    app_GPi_c1    = np.delete(app_GPi_c1, 0, 1)
    app_GPi_c2    = np.delete(app_GPi_c2, 0, 1)

    eps = np.arange(from_here, to_here)

    # D1-MSN mean ± std
    MSN_d1_c1     = np.mean(app_MSN_d1_c1, axis=1)
    MSN_d1_c1_err = np.std(app_MSN_d1_c1,  axis=1)
    MSN_d1_c2     = np.mean(app_MSN_d1_c2, axis=1)
    MSN_d1_c2_err = np.std(app_MSN_d1_c2,  axis=1)

    ax[0].plot(eps, MSN_d1_c1, '-',  color=my_color, label='Ch1 (mean)')
    ax[0].fill_between(eps, MSN_d1_c1 - MSN_d1_c1_err, MSN_d1_c1 + MSN_d1_c1_err,
                       alpha=0.2, color=my_color)
    ax[0].plot(eps, MSN_d1_c2, '-.', color=my_color, label='Ch2 (mean)')
    ax[0].fill_between(eps, MSN_d1_c2 - MSN_d1_c2_err, MSN_d1_c2 + MSN_d1_c2_err,
                       alpha=0.2, color=my_color)

    # GPi mean ± std
    GPi_c1     = np.mean(app_GPi_c1, axis=1)
    GPi_c1_err = np.std(app_GPi_c1,  axis=1)
    GPi_c2     = np.mean(app_GPi_c2, axis=1)
    GPi_c2_err = np.std(app_GPi_c2,  axis=1)

    ax[1].plot(eps, GPi_c1, '-',  color=my_color)
    ax[1].fill_between(eps, GPi_c1 - GPi_c1_err, GPi_c1 + GPi_c1_err,
                       alpha=0.2, color=my_color)
    ax[1].plot(eps, GPi_c2, '-.', color=my_color)
    ax[1].fill_between(eps, GPi_c2 - GPi_c2_err, GPi_c2 + GPi_c2_err,
                       alpha=0.2, color=my_color)

    ax[0].set_ylim([-0.1, 2.8])
    ax[0].set_title('D1-MSN\n\n$\\lambda=$' + my_lambda + ', $\\kappa$=' + my_kappa, fontsize=18)
    ax[1].set_title('GPi', fontsize=18)

    for a in ax:
        a.spines['right'].set_visible(False)
        a.spines['top'].set_visible(False)
        a.set_ylabel('Mean firing rate [Hz]', fontsize=16)
        a.set_xlabel('Episodes', fontsize=16)
        a.tick_params(axis='both', which='major', labelsize=14)

    ax[0].legend(loc='center left', fontsize=16)
    plt.subplots_adjust(hspace=0.35)

    if save:
        plt.savefig('fig_A_' + str(time.time()) + '.png', format='png', dpi=300)
    plt.show()


def Figure_4_D(save=True):
    """
    Figure 4D: CS+ conditioning — κ=2.0 vs κ=4.0 comparison.
    Generates 3 separate figures, one per λ value (0.0, 0.2, 0.4).
    """
    lambdas = [
        ('0.0', 'run4_k_2_l_002_m_1', 'run4_k_4_l_002_m_1'),
        ('0.2', 'run4_k_2_l_02_m_1',  'run4_k_4_l_02_m_1'),
        ('0.4', 'run4_k_2_l_04_m_1',  'run4_k_4_l_04_m_1'),
    ]
    my_color   = ['midnightblue', 'steelblue']   # κ=2, κ=4 (D1 / GPi Ch1)
    my_color_b = ['crimson',      'darkred']      # κ=2, κ=4 (D2 / GPi Ch2)
    up_to_ep = 100
    c1, c2 = 1, 2

    for lam, folder_k2, folder_k4 in lambdas:
        fig, axes = plt.subplots(2, figsize=(8, 10))
        ax_msn = axes[0]
        ax_gpi = axes[1]

        for k, folder in enumerate([folder_k2, folder_k4]):
            base = './' + folder + '/'
            app_MSN_d1_c1 = np.zeros((up_to_ep, 1))
            app_MSN_d1_c2 = np.zeros((up_to_ep, 1))
            app_MSN_d2_c1 = np.zeros((up_to_ep, 1))
            app_MSN_d2_c2 = np.zeros((up_to_ep, 1))
            app_GPi_c1    = np.zeros((up_to_ep, 1))
            app_GPi_c2    = np.zeros((up_to_ep, 1))
            my_kappa = ''

            for m in np.arange(10):
                my_path = base + 'log' + str(m) + '/'
                try:
                    MSN_d1_fr = np.loadtxt(my_path + 'MSN_d1_fr.txt')
                    MSN_d2_fr = np.loadtxt(my_path + 'MSN_d2_fr.txt')
                    GPi_fr    = np.loadtxt(my_path + 'GPi_fr.txt')

                    app_MSN_d1_c1 = np.append(app_MSN_d1_c1, MSN_d1_fr[:up_to_ep, c1].reshape(up_to_ep, 1), axis=1)
                    app_MSN_d1_c2 = np.append(app_MSN_d1_c2, MSN_d1_fr[:up_to_ep, c2].reshape(up_to_ep, 1), axis=1)
                    app_MSN_d2_c1 = np.append(app_MSN_d2_c1, MSN_d2_fr[:up_to_ep, c1].reshape(up_to_ep, 1), axis=1)
                    app_MSN_d2_c2 = np.append(app_MSN_d2_c2, MSN_d2_fr[:up_to_ep, c2].reshape(up_to_ep, 1), axis=1)
                    app_GPi_c1    = np.append(app_GPi_c1,    GPi_fr[:up_to_ep, c1].reshape(up_to_ep, 1),    axis=1)
                    app_GPi_c2    = np.append(app_GPi_c2,    GPi_fr[:up_to_ep, c2].reshape(up_to_ep, 1),    axis=1)

                    params = np.loadtxt(my_path + 'plot_labels.txt', dtype=str)
                    my_kappa = params[1].split(',')[1]
                except Exception:
                    pass

            # Remove initialisation column
            app_MSN_d1_c1 = np.delete(app_MSN_d1_c1, 0, 1)
            app_MSN_d1_c2 = np.delete(app_MSN_d1_c2, 0, 1)
            app_MSN_d2_c1 = np.delete(app_MSN_d2_c1, 0, 1)
            app_MSN_d2_c2 = np.delete(app_MSN_d2_c2, 0, 1)
            app_GPi_c1    = np.delete(app_GPi_c1, 0, 1)
            app_GPi_c2    = np.delete(app_GPi_c2, 0, 1)

            eps = np.arange(1, up_to_ep + 1)

            # D1-MSN (blue shades)
            d1_c1 = np.mean(app_MSN_d1_c1, axis=1); d1_c1_err = np.std(app_MSN_d1_c1, axis=1)
            d1_c2 = np.mean(app_MSN_d1_c2, axis=1); d1_c2_err = np.std(app_MSN_d1_c2, axis=1)
            ax_msn.plot(eps, d1_c1, '-',  color=my_color[k],   label='$\\kappa$=' + my_kappa + ', D1-MSN')
            ax_msn.fill_between(eps, d1_c1 - d1_c1_err, d1_c1 + d1_c1_err, alpha=0.15, color=my_color[k])
            ax_msn.plot(eps, d1_c2, '-.', color=my_color[k])
            ax_msn.fill_between(eps, d1_c2 - d1_c2_err, d1_c2 + d1_c2_err, alpha=0.15, color=my_color[k])

            # D2-MSN (red shades)
            d2_c1 = np.mean(app_MSN_d2_c1, axis=1); d2_c1_err = np.std(app_MSN_d2_c1, axis=1)
            d2_c2 = np.mean(app_MSN_d2_c2, axis=1); d2_c2_err = np.std(app_MSN_d2_c2, axis=1)
            ax_msn.plot(eps, d2_c1, '-',  color=my_color_b[k], label='$\\kappa$=' + my_kappa + ', D2-MSN')
            ax_msn.fill_between(eps, d2_c1 - d2_c1_err, d2_c1 + d2_c1_err, alpha=0.15, color=my_color_b[k])
            ax_msn.plot(eps, d2_c2, '-.', color=my_color_b[k])
            ax_msn.fill_between(eps, d2_c2 - d2_c2_err, d2_c2 + d2_c2_err, alpha=0.15, color=my_color_b[k])

            # GPi Ch1 (blue) / Ch2 (red)
            g_c1 = np.mean(app_GPi_c1, axis=1); g_c1_err = np.std(app_GPi_c1, axis=1)
            g_c2 = np.mean(app_GPi_c2, axis=1); g_c2_err = np.std(app_GPi_c2, axis=1)
            ax_gpi.plot(eps, g_c1, '-',  color=my_color[k],   label='$\\kappa$=' + my_kappa)
            ax_gpi.fill_between(eps, g_c1 - g_c1_err, g_c1 + g_c1_err, alpha=0.15, color=my_color[k])
            ax_gpi.plot(eps, g_c2, '-.',  color=my_color_b[k], label='$\\kappa$=' + my_kappa)
            ax_gpi.fill_between(eps, g_c2 - g_c2_err, g_c2 + g_c2_err, alpha=0.15, color=my_color_b[k])

        # Legend for MSN panel
        legend_msn = [
            Line2D([0], [0], ls='-',  color='k',             label='Ch1 (mean)'),
            Line2D([0], [0], ls='-.', color='k',             label='Ch2 (mean)'),
            Line2D([0], [0], ls='',   marker='s', color=my_color[0],   markersize=7, label='$\\kappa$=2.0, D1-MSN'),
            Line2D([0], [0], ls='',   marker='s', color=my_color[1],   markersize=7, label='$\\kappa$=4.0, D1-MSN'),
            Line2D([0], [0], ls='',   marker='s', color=my_color_b[0], markersize=7, label='$\\kappa$=2.0, D2-MSN'),
            Line2D([0], [0], ls='',   marker='s', color=my_color_b[1], markersize=7, label='$\\kappa$=4.0, D2-MSN'),
        ]
        leg = ax_msn.legend(handles=legend_msn, loc=(0.05, 0.58), fontsize=16, frameon=False, labelspacing=0.02)
        leg._legend_box.align = 'left'

        # Legend for GPi panel
        legend_gpi = [
            Line2D([0], [0], ls='-',  color='k',             label='Ch1 (mean)'),
            Line2D([0], [0], ls='',   marker='s', color=my_color[0],   markersize=7, label='$\\kappa$=2.0'),
            Line2D([0], [0], ls='',   marker='s', color=my_color[1],   markersize=7, label='$\\kappa$=4.0'),
            Line2D([0], [0], ls='-.', color='k',             label='Ch2 (mean)'),
            Line2D([0], [0], ls='',   marker='s', color=my_color_b[0], markersize=7, label='$\\kappa$=2.0'),
            Line2D([0], [0], ls='',   marker='s', color=my_color_b[1], markersize=7, label='$\\kappa$=4.0'),
        ]
        leg2 = ax_gpi.legend(handles=legend_gpi, loc=(0.05, 0.12), fontsize=16, frameon=False, labelspacing=0.02)
        leg2._legend_box.align = 'left'

        ax_msn.set_ylim([-0.5, 2.0])
        ax_gpi.set_ylim([42, 87])
        ax_msn.set_title('MSN (CS+)\n\n$\\lambda=$' + lam, fontsize=18)
        ax_gpi.set_title('GPi', fontsize=18)

        for ax in [ax_msn, ax_gpi]:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_xlabel('Episodes', fontsize=16)
            ax.set_ylabel('Mean firing rate [Hz]', fontsize=16)
            ax.tick_params(axis='both', which='major', labelsize=14)

        plt.subplots_adjust(hspace=0.35)
        if save:
            plt.savefig('fig_4D_lambda' + lam + '_' + str(time.time()) + '.png', format='png', dpi=300)
        plt.show()


def Figure_5_A(save=True):
    """
    Figure 5A: CS- discrimination — κ=2.0 vs κ=4.0 for λ=0.0, 0.1, 0.2, 0.3.
    Generates 4 separate figures (one per λ), episodes 131–225.
    D2-MSN in red shades (Ch2 stimulated), D1-MSN in blue shades.
    """
    lambdas = [
        ('0.0', 'run4_k_2_l_002_m_1', 'run4_k_4_l_002_m_1'),
        ('0.1', 'run4_k_2_l_01_m_1',  'run4_k_4_l_01_m_1'),
        ('0.2', 'run4_k_2_l_02_m_1',  'run4_k_4_l_02_m_1'),
        ('0.3', 'run4_k_2_l_03_m_1',  'run4_k_4_l_03_m_1'),
    ]
    my_color   = ['crimson',      'darkred']      # κ=2, κ=4 (D2-MSN / GPi Ch2)
    my_color_b = ['midnightblue', 'steelblue']    # κ=2, κ=4 (D1-MSN / GPi Ch1)
    from_here = 131
    to_here   = 225
    c1, c2 = 1, 2

    for lam, folder_k2, folder_k4 in lambdas:
        fig, ax = plt.subplots(2, figsize=(8, 10))

        for k, folder in enumerate([folder_k2, folder_k4]):
            base = './' + folder + '/'
            n = to_here - from_here
            app_MSN_d1_c1 = np.zeros((n, 1))
            app_MSN_d1_c2 = np.zeros((n, 1))
            app_MSN_d2_c1 = np.zeros((n, 1))
            app_MSN_d2_c2 = np.zeros((n, 1))
            app_GPi_c1    = np.zeros((n, 1))
            app_GPi_c2    = np.zeros((n, 1))
            my_kappa = ''

            for m in np.arange(20):
                my_path = base + 'log' + str(m) + '/'
                try:
                    MSN_d1_fr = np.loadtxt(my_path + 'MSN_d1_fr.txt')
                    MSN_d2_fr = np.loadtxt(my_path + 'MSN_d2_fr.txt')
                    GPi_fr    = np.loadtxt(my_path + 'GPi_fr.txt')

                    app_MSN_d1_c1 = np.append(app_MSN_d1_c1, MSN_d1_fr[from_here:to_here, c1].reshape(n, 1), axis=1)
                    app_MSN_d1_c2 = np.append(app_MSN_d1_c2, MSN_d1_fr[from_here:to_here, c2].reshape(n, 1), axis=1)
                    app_MSN_d2_c1 = np.append(app_MSN_d2_c1, MSN_d2_fr[from_here:to_here, c1].reshape(n, 1), axis=1)
                    app_MSN_d2_c2 = np.append(app_MSN_d2_c2, MSN_d2_fr[from_here:to_here, c2].reshape(n, 1), axis=1)
                    app_GPi_c1    = np.append(app_GPi_c1,    GPi_fr[from_here:to_here, c1].reshape(n, 1),    axis=1)
                    app_GPi_c2    = np.append(app_GPi_c2,    GPi_fr[from_here:to_here, c2].reshape(n, 1),    axis=1)

                    params = np.loadtxt(my_path + 'plot_labels.txt', dtype=str)
                    my_kappa = params[1].split(',')[1]
                except Exception:
                    pass

            for arr in [app_MSN_d1_c1, app_MSN_d1_c2, app_MSN_d2_c1,
                        app_MSN_d2_c2, app_GPi_c1, app_GPi_c2]:
                arr = np.delete(arr, 0, 1)

            app_MSN_d1_c1 = np.delete(app_MSN_d1_c1, 0, 1)
            app_MSN_d1_c2 = np.delete(app_MSN_d1_c2, 0, 1)
            app_MSN_d2_c1 = np.delete(app_MSN_d2_c1, 0, 1)
            app_MSN_d2_c2 = np.delete(app_MSN_d2_c2, 0, 1)
            app_GPi_c1    = np.delete(app_GPi_c1, 0, 1)
            app_GPi_c2    = np.delete(app_GPi_c2, 0, 1)

            eps = np.arange(from_here, to_here)

            d1_c1 = np.mean(app_MSN_d1_c1, axis=1); d1_c1_err = np.std(app_MSN_d1_c1, axis=1)
            d1_c2 = np.mean(app_MSN_d1_c2, axis=1); d1_c2_err = np.std(app_MSN_d1_c2, axis=1)
            d2_c1 = np.mean(app_MSN_d2_c1, axis=1); d2_c1_err = np.std(app_MSN_d2_c1, axis=1)
            d2_c2 = np.mean(app_MSN_d2_c2, axis=1); d2_c2_err = np.std(app_MSN_d2_c2, axis=1)

            # D2-MSN (red) — Ch2 stimulated
            ax[0].plot(eps, d2_c2, '-.', color=my_color[k])
            ax[0].fill_between(eps, d2_c2 - d2_c2_err, d2_c2 + d2_c2_err, alpha=0.2, color=my_color[k])
            ax[0].plot(eps, d2_c1, '-',  color=my_color[k])
            ax[0].fill_between(eps, d2_c1 - d2_c1_err, d2_c1 + d2_c1_err, alpha=0.2, color=my_color[k])

            # D1-MSN (blue)
            ax[0].plot(eps, d1_c2, '-.', color=my_color_b[k])
            ax[0].fill_between(eps, d1_c2 - d1_c2_err, d1_c2 + d1_c2_err, alpha=0.2, color=my_color_b[k])
            ax[0].plot(eps, d1_c1, '-',  color=my_color_b[k])
            ax[0].fill_between(eps, d1_c1 - d1_c1_err, d1_c1 + d1_c1_err, alpha=0.2, color=my_color_b[k])

            g_c1 = np.mean(app_GPi_c1, axis=1); g_c1_err = np.std(app_GPi_c1, axis=1)
            g_c2 = np.mean(app_GPi_c2, axis=1); g_c2_err = np.std(app_GPi_c2, axis=1)

            ax[1].plot(eps, g_c1, '-',  color=my_color_b[k])
            ax[1].fill_between(eps, g_c1 - g_c1_err, g_c1 + g_c1_err, alpha=0.2, color=my_color_b[k])
            ax[1].plot(eps, g_c2, '-.', color=my_color[k])
            ax[1].fill_between(eps, g_c2 - g_c2_err, g_c2 + g_c2_err, alpha=0.2, color=my_color[k])

        legend_msn = [
            Line2D([0], [0], ls='-',  color='k',              label='Ch1 (mean)'),
            Line2D([0], [0], ls='-.', color='k',              label='Ch2 (mean)'),
            Line2D([0], [0], ls='',   marker='s', color=my_color_b[0], markersize=7, label='$\\kappa$=2.0, D1-MSN'),
            Line2D([0], [0], ls='',   marker='s', color=my_color_b[1], markersize=7, label='$\\kappa$=4.0, D1-MSN'),
            Line2D([0], [0], ls='',   marker='s', color=my_color[0],   markersize=7, label='$\\kappa$=2.0, D2-MSN'),
            Line2D([0], [0], ls='',   marker='s', color=my_color[1],   markersize=7, label='$\\kappa$=4.0, D2-MSN'),
        ]
        leg0 = ax[0].legend(handles=legend_msn, loc=(0.05, 0.58), fontsize=16, frameon=False, labelspacing=0.02)
        leg0._legend_box.align = 'left'

        legend_gpi = [
            Line2D([0], [0], ls='-',  color='k',              label='Ch1 (mean)'),
            Line2D([0], [0], ls='',   marker='s', color=my_color_b[0], markersize=7, label='$\\kappa$=2.0'),
            Line2D([0], [0], ls='',   marker='s', color=my_color_b[1], markersize=7, label='$\\kappa$=4.0'),
            Line2D([0], [0], ls='-.', color='k',              label='Ch2 (mean)'),
            Line2D([0], [0], ls='',   marker='s', color=my_color[0],   markersize=7, label='$\\kappa$=2.0'),
            Line2D([0], [0], ls='',   marker='s', color=my_color[1],   markersize=7, label='$\\kappa$=4.0'),
        ]
        leg1 = ax[1].legend(handles=legend_gpi, loc=(0.16, 0.001), fontsize=16, frameon=False,
                            labelspacing=0.02, ncol=2)
        leg1._legend_box.align = 'left'

        ax[0].set_ylim([-0.5, 10])
        ax[1].set_ylim([50, 85])
        ax[0].set_title('MSN (CS-)\n\n$\\lambda=$' + lam, fontsize=18)
        ax[1].set_title('GPi', fontsize=18)

        for a in ax:
            a.spines['right'].set_visible(False)
            a.spines['top'].set_visible(False)
            a.set_ylabel('Mean firing rate [Hz]', fontsize=16)
            a.set_xlabel('Episodes', fontsize=16)
            a.tick_params(axis='both', which='major', labelsize=14)

        plt.subplots_adjust(hspace=0.35)
        if save:
            plt.savefig('fig_5A_lambda' + lam + '_' + str(time.time()) + '.png', format='png', dpi=300)
        plt.show()


def Figure_6_AB(save=True):
    """
    Figure 6 A&B: CS- zoomed MSN — same data as Figure_5_A but single MSN panel
    with tight ylim [-0.5, 1.5] Hz. 4 separate files, one per λ (0.0, 0.1, 0.2, 0.3).
    """
    lambdas = [
        ('0.0', 'run4_k_2_l_002_m_1', 'run4_k_4_l_002_m_1'),
        ('0.1', 'run4_k_2_l_01_m_1',  'run4_k_4_l_01_m_1'),
        ('0.2', 'run4_k_2_l_02_m_1',  'run4_k_4_l_02_m_1'),
        ('0.3', 'run4_k_2_l_03_m_1',  'run4_k_4_l_03_m_1'),
    ]
    my_color   = ['crimson',      'darkred']      # κ=2, κ=4 (D2-MSN)
    my_color_b = ['midnightblue', 'steelblue']    # κ=2, κ=4 (D1-MSN)
    from_here = 131
    to_here   = 225
    c1, c2 = 1, 2

    for lam, folder_k2, folder_k4 in lambdas:
        fig, ax = plt.subplots(1, figsize=(8, 5))

        for k, folder in enumerate([folder_k2, folder_k4]):
            base = './' + folder + '/'
            n = to_here - from_here
            app_MSN_d1_c1 = np.zeros((n, 1))
            app_MSN_d1_c2 = np.zeros((n, 1))
            app_MSN_d2_c1 = np.zeros((n, 1))
            app_MSN_d2_c2 = np.zeros((n, 1))
            my_kappa = ''

            for m in np.arange(20):
                my_path = base + 'log' + str(m) + '/'
                try:
                    MSN_d1_fr = np.loadtxt(my_path + 'MSN_d1_fr.txt')
                    MSN_d2_fr = np.loadtxt(my_path + 'MSN_d2_fr.txt')

                    app_MSN_d1_c1 = np.append(app_MSN_d1_c1, MSN_d1_fr[from_here:to_here, c1].reshape(n, 1), axis=1)
                    app_MSN_d1_c2 = np.append(app_MSN_d1_c2, MSN_d1_fr[from_here:to_here, c2].reshape(n, 1), axis=1)
                    app_MSN_d2_c1 = np.append(app_MSN_d2_c1, MSN_d2_fr[from_here:to_here, c1].reshape(n, 1), axis=1)
                    app_MSN_d2_c2 = np.append(app_MSN_d2_c2, MSN_d2_fr[from_here:to_here, c2].reshape(n, 1), axis=1)

                    params = np.loadtxt(my_path + 'plot_labels.txt', dtype=str)
                    my_kappa = params[1].split(',')[1]
                except Exception:
                    pass

            app_MSN_d1_c1 = np.delete(app_MSN_d1_c1, 0, 1)
            app_MSN_d1_c2 = np.delete(app_MSN_d1_c2, 0, 1)
            app_MSN_d2_c1 = np.delete(app_MSN_d2_c1, 0, 1)
            app_MSN_d2_c2 = np.delete(app_MSN_d2_c2, 0, 1)

            eps = np.arange(from_here, to_here)

            d1_c1 = np.mean(app_MSN_d1_c1, axis=1); d1_c1_err = np.std(app_MSN_d1_c1, axis=1)
            d1_c2 = np.mean(app_MSN_d1_c2, axis=1); d1_c2_err = np.std(app_MSN_d1_c2, axis=1)
            d2_c1 = np.mean(app_MSN_d2_c1, axis=1); d2_c1_err = np.std(app_MSN_d2_c1, axis=1)
            d2_c2 = np.mean(app_MSN_d2_c2, axis=1); d2_c2_err = np.std(app_MSN_d2_c2, axis=1)

            # D2-MSN (red)
            ax.plot(eps, d2_c2, '-.', color=my_color[k])
            ax.fill_between(eps, d2_c2 - d2_c2_err, d2_c2 + d2_c2_err, alpha=0.2, color=my_color[k])
            ax.plot(eps, d2_c1, '-',  color=my_color[k])
            ax.fill_between(eps, d2_c1 - d2_c1_err, d2_c1 + d2_c1_err, alpha=0.2, color=my_color[k])

            # D1-MSN (blue)
            ax.plot(eps, d1_c2, '-.', color=my_color_b[k])
            ax.fill_between(eps, d1_c2 - d1_c2_err, d1_c2 + d1_c2_err, alpha=0.2, color=my_color_b[k])
            ax.plot(eps, d1_c1, '-',  color=my_color_b[k])
            ax.fill_between(eps, d1_c1 - d1_c1_err, d1_c1 + d1_c1_err, alpha=0.2, color=my_color_b[k])

        legend_msn = [
            Line2D([0], [0], ls='-',  color='k',              label='Ch1 (mean)'),
            Line2D([0], [0], ls='-.', color='k',              label='Ch2 (mean)'),
            Line2D([0], [0], ls='',   marker='s', color=my_color_b[0], markersize=7, label='$\\kappa$=2.0, D1-MSN'),
            Line2D([0], [0], ls='',   marker='s', color=my_color_b[1], markersize=7, label='$\\kappa$=4.0, D1-MSN'),
            Line2D([0], [0], ls='',   marker='s', color=my_color[0],   markersize=7, label='$\\kappa$=2.0, D2-MSN'),
            Line2D([0], [0], ls='',   marker='s', color=my_color[1],   markersize=7, label='$\\kappa$=4.0, D2-MSN'),
        ]
        leg = ax.legend(handles=legend_msn, loc='lower left', fontsize=16, frameon=False,
                        labelspacing=0.02, ncol=2)
        leg._legend_box.align = 'left'

        ax.set_ylim([-0.5, 1.5])
        ax.set_title('MSN (CS-)\n\n$\\lambda=$' + lam, fontsize=18)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_ylabel('Mean firing rate [Hz]', fontsize=16)
        ax.set_xlabel('Episodes', fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=14)

        plt.tight_layout()
        if save:
            plt.savefig('fig_6AB_lambda' + lam + '_' + str(time.time()) + '.png', format='png', dpi=300)
        plt.show()


def Figure_6_C(save=True):
    """
    Figure 6C: Full learning trajectory (CS+ + GT + CS-), κ=4.0, λ=0.2/0.3/0.4.
    Single figure with 2 panels: MSN (top) and GPi (bottom) with inset zoom.
    D1-MSN Ch1 (dots) and D2-MSN Ch2 (plus) per λ, episodes 0–225.
    """
    folders = ['run4_k_4_l_02_m_1', 'run4_k_4_l_03_m_1', 'run4_k_4_l_04_m_1']
    lambdas_label = ['0.2', '0.3', '0.4']
    my_color = ['royalblue', 'crimson', 'mediumorchid']  # one per λ
    from_here = 0
    to_here   = 225
    c1, c2 = 1, 2

    fig, ax = plt.subplots(2, figsize=(12, 12))
    axin = ax[1].inset_axes([0.7, 0.1, 0.35, 0.25])

    my_kappa = '4.0'

    for k, (folder, lam) in enumerate(zip(folders, lambdas_label)):
        base = './' + folder + '/'
        n = to_here - from_here
        app_MSN_d1_c1 = np.zeros((n, 1))
        app_MSN_d2_c2 = np.zeros((n, 1))
        app_GPi_c1    = np.zeros((n, 1))
        app_GPi_c2    = np.zeros((n, 1))

        for m in np.arange(20):
            my_path = base + 'log' + str(m) + '/'
            try:
                MSN_d1_fr = np.loadtxt(my_path + 'MSN_d1_fr.txt')
                MSN_d2_fr = np.loadtxt(my_path + 'MSN_d2_fr.txt')
                GPi_fr    = np.loadtxt(my_path + 'GPi_fr.txt')

                app_MSN_d1_c1 = np.append(app_MSN_d1_c1, MSN_d1_fr[from_here:to_here, c1].reshape(n, 1), axis=1)
                app_MSN_d2_c2 = np.append(app_MSN_d2_c2, MSN_d2_fr[from_here:to_here, c2].reshape(n, 1), axis=1)
                app_GPi_c1    = np.append(app_GPi_c1,    GPi_fr[from_here:to_here, c1].reshape(n, 1),    axis=1)
                app_GPi_c2    = np.append(app_GPi_c2,    GPi_fr[from_here:to_here, c2].reshape(n, 1),    axis=1)
            except Exception:
                pass

        app_MSN_d1_c1 = np.delete(app_MSN_d1_c1, 0, 1)
        app_MSN_d2_c2 = np.delete(app_MSN_d2_c2, 0, 1)
        app_GPi_c1    = np.delete(app_GPi_c1, 0, 1)
        app_GPi_c2    = np.delete(app_GPi_c2, 0, 1)

        eps = np.arange(from_here, to_here)

        d1_c1     = np.mean(app_MSN_d1_c1, axis=1); d1_c1_err = np.std(app_MSN_d1_c1, axis=1)
        d2_c2     = np.mean(app_MSN_d2_c2, axis=1); d2_c2_err = np.std(app_MSN_d2_c2, axis=1)
        g_c1      = np.mean(app_GPi_c1,    axis=1); g_c1_err  = np.std(app_GPi_c1,    axis=1)
        g_c2      = np.mean(app_GPi_c2,    axis=1); g_c2_err  = np.std(app_GPi_c2,    axis=1)

        col = my_color[k]
        ax[0].plot(eps, d1_c1, '.', color=col, label='D1-MSN (Ch1), $\\lambda$=' + lam, alpha=0.6, markeredgewidth=1.1)
        ax[0].fill_between(eps, d1_c1 - d1_c1_err, d1_c1 + d1_c1_err, alpha=0.1, color=col)
        ax[0].plot(eps, d2_c2, '+', color=col, label='D2-MSN (Ch2), $\\lambda$=' + lam, alpha=0.6, markeredgewidth=1.1)
        ax[0].fill_between(eps, d2_c2 - d2_c2_err, d2_c2 + d2_c2_err, alpha=0.1, color=col)

        ax[1].plot(eps, g_c1, '.', color=col, label='GPi (Ch1), $\\lambda$=' + lam, alpha=0.6, markeredgewidth=1.1)
        ax[1].fill_between(eps, g_c1 - g_c1_err, g_c1 + g_c1_err, alpha=0.1, color=col)
        ax[1].plot(eps, g_c2, '+', color=col, label='GPi (Ch2), $\\lambda$=' + lam, alpha=0.6, markeredgewidth=1.1)
        ax[1].fill_between(eps, g_c2 - g_c2_err, g_c2 + g_c2_err, alpha=0.1, color=col)

        axin.plot(eps, g_c1, '.', color=col, alpha=0.6, markeredgewidth=1.1)
        axin.plot(eps, g_c2, '+', color=col, alpha=0.6, markeredgewidth=1.1)

    # Phase boundary lines
    for a in ax:
        ylims = [(-1, 7.5), (40, 85.5)]
    ax[0].plot([100, 100], [-1,  7.5], 'k-.', alpha=0.5)
    ax[0].plot([130, 130], [-1,  7.5], 'k-.', alpha=0.5)
    ax[0].plot([225, 225], [-1,  7.5], 'k-.', alpha=0.5)
    ax[1].plot([100, 100], [40, 85.5], 'k-.', alpha=0.5)
    ax[1].plot([130, 130], [40, 85.5], 'k-.', alpha=0.5)
    ax[1].plot([225, 225], [40, 85.5], 'k-.', alpha=0.5)


    ax[0].legend(loc='upper left', fontsize=16)
    ax[1].legend(loc='lower left', fontsize=16)

    ax[0].set_title('MSN\n\n$\\kappa$=' + my_kappa, fontsize=18)
    ax[1].set_title('GPi', fontsize=18)

    for a in ax:
        a.spines['right'].set_visible(False)
        a.spines['top'].set_visible(False)
        a.set_ylabel('Mean firing rate [Hz]', fontsize=16)
        a.tick_params(axis='both', which='major', labelsize=14)
    ax[1].set_xlabel('Episodes', fontsize=16)

    axin.set_xlim(153, 220)
    axin.set_ylim(60, 77)
    ax[1].indicate_inset_zoom(axin)

    plt.subplots_adjust(hspace=0.35)
    if save:
        plt.savefig('fig_6C_' + str(time.time()) + '.png', format='png', dpi=300)
    plt.show()


def Figure_7_D(save=True):
    """
    Figure 7D: Full learning trajectory (CS+ + GT + CS-), λ=0.2, κ=4.0, η=1.0/2.0/4.0.
    Same structure as Figure_6_C but varying modulation (η) instead of overlap (λ).
    """
    folders    = ['run4_k_4_l_02_m_1', 'run4_k_4_l_02_m_2', 'run4_k_4_l_02_m_4']
    eta_labels = ['1.0', '2.0', '4.0']
    my_color   = ['royalblue', 'crimson', 'mediumorchid']
    from_here = 0
    to_here   = 225
    c1, c2 = 1, 2

    fig, ax = plt.subplots(2, figsize=(12, 12))
    axin = ax[1].inset_axes([0.7, 0.1, 0.35, 0.25])

    for k, (folder, eta) in enumerate(zip(folders, eta_labels)):
        base = './' + folder + '/'
        n = to_here - from_here
        app_MSN_d1_c1 = np.zeros((n, 1))
        app_MSN_d2_c2 = np.zeros((n, 1))
        app_GPi_c1    = np.zeros((n, 1))
        app_GPi_c2    = np.zeros((n, 1))

        for m in np.arange(20):
            my_path = base + 'log' + str(m) + '/'
            try:
                MSN_d1_fr = np.loadtxt(my_path + 'MSN_d1_fr.txt')
                MSN_d2_fr = np.loadtxt(my_path + 'MSN_d2_fr.txt')
                GPi_fr    = np.loadtxt(my_path + 'GPi_fr.txt')

                app_MSN_d1_c1 = np.append(app_MSN_d1_c1, MSN_d1_fr[from_here:to_here, c1].reshape(n, 1), axis=1)
                app_MSN_d2_c2 = np.append(app_MSN_d2_c2, MSN_d2_fr[from_here:to_here, c2].reshape(n, 1), axis=1)
                app_GPi_c1    = np.append(app_GPi_c1,    GPi_fr[from_here:to_here, c1].reshape(n, 1),    axis=1)
                app_GPi_c2    = np.append(app_GPi_c2,    GPi_fr[from_here:to_here, c2].reshape(n, 1),    axis=1)
            except Exception:
                pass

        app_MSN_d1_c1 = np.delete(app_MSN_d1_c1, 0, 1)
        app_MSN_d2_c2 = np.delete(app_MSN_d2_c2, 0, 1)
        app_GPi_c1    = np.delete(app_GPi_c1, 0, 1)
        app_GPi_c2    = np.delete(app_GPi_c2, 0, 1)

        eps = np.arange(from_here, to_here)
        col = my_color[k]

        d1_c1 = np.mean(app_MSN_d1_c1, axis=1); d1_c1_err = np.std(app_MSN_d1_c1, axis=1)
        d2_c2 = np.mean(app_MSN_d2_c2, axis=1); d2_c2_err = np.std(app_MSN_d2_c2, axis=1)
        g_c1  = np.mean(app_GPi_c1,    axis=1); g_c1_err  = np.std(app_GPi_c1,    axis=1)
        g_c2  = np.mean(app_GPi_c2,    axis=1); g_c2_err  = np.std(app_GPi_c2,    axis=1)

        ax[0].plot(eps, d1_c1, '.', color=col, label='D1-MSN (Ch1), $\\eta$=' + eta, alpha=0.6, markeredgewidth=1.1)
        ax[0].fill_between(eps, d1_c1 - d1_c1_err, d1_c1 + d1_c1_err, alpha=0.1, color=col)
        ax[0].plot(eps, d2_c2, '+', color=col, label='D2-MSN (Ch2), $\\eta$=' + eta, alpha=0.6, markeredgewidth=1.1)
        ax[0].fill_between(eps, d2_c2 - d2_c2_err, d2_c2 + d2_c2_err, alpha=0.1, color=col)

        ax[1].plot(eps, g_c1, '.', color=col, label='GPi (Ch1), $\\eta$=' + eta, alpha=0.6, markeredgewidth=1.1)
        ax[1].fill_between(eps, g_c1 - g_c1_err, g_c1 + g_c1_err, alpha=0.1, color=col)
        ax[1].plot(eps, g_c2, '+', color=col, label='GPi (Ch2), $\\eta$=' + eta, alpha=0.6, markeredgewidth=1.1)
        ax[1].fill_between(eps, g_c2 - g_c2_err, g_c2 + g_c2_err, alpha=0.1, color=col)

        axin.plot(eps, g_c1, '.', color=col, alpha=0.6, markeredgewidth=1.1)
        axin.plot(eps, g_c2, '+', color=col, alpha=0.6, markeredgewidth=1.1)

    # Phase boundary lines
    ax[0].plot([100, 100], [-1,  7.5], 'k-.', alpha=0.5)
    ax[0].plot([130, 130], [-1,  7.5], 'k-.', alpha=0.5)
    ax[0].plot([225, 225], [-1,  7.5], 'k-.', alpha=0.5)
    ax[1].plot([100, 100], [40, 85.5], 'k-.', alpha=0.5)
    ax[1].plot([130, 130], [40, 85.5], 'k-.', alpha=0.5)
    ax[1].plot([225, 225], [40, 85.5], 'k-.', alpha=0.5)

    ax[0].legend(loc='upper left', fontsize=16)
    ax[1].legend(loc='lower left', fontsize=16)

    ax[0].set_title('MSN\n\n$\\lambda$=0.2, $\\kappa$=4.0', fontsize=18)
    ax[1].set_title('GPi', fontsize=18)

    for a in ax:
        a.spines['right'].set_visible(False)
        a.spines['top'].set_visible(False)
        a.set_ylabel('Mean firing rate [Hz]', fontsize=16)
        a.tick_params(axis='both', which='major', labelsize=14)
    ax[1].set_xlabel('Episodes', fontsize=16)

    axin.set_xlim(153, 220)
    axin.set_ylim(60, 77)
    ax[1].indicate_inset_zoom(axin)

    plt.subplots_adjust(hspace=0.35)
    if save:
        plt.savefig('fig_7D_' + str(time.time()) + '.png', format='png', dpi=300)
    plt.show()


def _plot_individual_trajectory(folder, epis, logid, gen_interv, msize, title, save, fname):
    """
    Helper for Figure 8 A/B/D: full trajectory from a single log, open-circle markers,
    generalization test markers, phase boundary lines.
    """
    my_color   = 'crimson'   # D2-MSN / GPi Ch2
    my_color_b = 'steelblue' # D1-MSN / GPi Ch1
    from_here, to_here = 0, epis
    c1, c2 = 1, 2
    n = to_here - from_here

    my_path = './' + folder + '/log' + str(logid) + '/'
    MSN_d1_fr = np.loadtxt(my_path + 'MSN_d1_fr.txt')
    MSN_d2_fr = np.loadtxt(my_path + 'MSN_d2_fr.txt')
    GPi_fr    = np.loadtxt(my_path + 'GPi_fr.txt')
    params    = np.loadtxt(my_path + 'plot_labels.txt', dtype=str)
    my_lambda = params[0].split(',')[1]
    my_kappa  = params[1].split(',')[1]

    # Clamp to actual data length in case simulation ended early
    to_here = min(to_here, MSN_d1_fr.shape[0])
    eps = np.arange(from_here, to_here)

    d1_c1 = MSN_d1_fr[from_here:to_here, c1]
    d2_c1 = MSN_d2_fr[from_here:to_here, c1]
    d1_c2 = MSN_d1_fr[from_here:to_here, c2]
    d2_c2 = MSN_d2_fr[from_here:to_here, c2]
    g_c1  = GPi_fr[from_here:to_here, c1]
    g_c2  = GPi_fr[from_here:to_here, c2]

    fig, ax = plt.subplots(2, figsize=(12, 12))

    kw = dict(alpha=0.7, markeredgewidth=1.1, markerfacecolor='None', markersize=msize)
    ax[0].plot(eps, d1_c1, 'o', markeredgecolor=my_color_b, label='D1-MSN (Ch1), $\\lambda$=' + my_lambda, **kw)
    ax[0].plot(eps, d2_c1, 'o', markeredgecolor=my_color,   label='D2-MSN (Ch1), $\\lambda$=' + my_lambda, **kw)
    ax[0].plot(eps, d2_c2, '+', markeredgecolor=my_color,   label='D2-MSN (Ch2), $\\lambda$=' + my_lambda, **kw)
    ax[0].plot(eps, d1_c2, '+', markeredgecolor=my_color_b, label='D1-MSN (Ch2), $\\lambda$=' + my_lambda, **kw)

    ax[1].plot(eps, g_c1, 'o', markeredgecolor=my_color_b, label='GPi (Ch1)',  **kw)
    ax[1].plot(eps, g_c2, '+', markeredgecolor=my_color,   label='GPi (Ch2)',  **kw)

    # Generalization test markers
    g0 = gen_interv[0]
    ax[0].plot([g0+5, g0+15], [d1_c1[g0+5], d1_c1[g0+15]], '*', color='k', label='Gen.Test (stim. A)', ms=msize, alpha=0.5)
    ax[0].plot([g0+10,g0+20], [d1_c1[g0+10],d1_c1[g0+20]], 's', color='k', label='Gen.Test (stim. B)', ms=msize, alpha=0.5)
    ax[0].plot([g0+5, g0+15], [d2_c2[g0+5], d2_c2[g0+15]], '*', color='k', ms=msize, alpha=0.5)
    ax[0].plot([g0+10,g0+20], [d2_c2[g0+10],d2_c2[g0+20]], 's', color='k', ms=msize, alpha=0.5)
    ax[1].plot([g0+5, g0+15], [g_c2[g0+5], g_c2[g0+15]], '*', color='k', label='Gen.Test (stim. A)', ms=msize, alpha=0.5)
    ax[1].plot([g0+10,g0+20], [g_c2[g0+10],g_c2[g0+20]], 's', color='k', label='Gen.Test (stim. B)', ms=msize, alpha=0.5)
    ax[1].plot([g0+5, g0+15], [g_c1[g0+5], g_c1[g0+15]], '*', color='k', ms=msize, alpha=0.5)
    ax[1].plot([g0+10,g0+20], [g_c1[g0+10],g_c1[g0+20]], 's', color='k', ms=msize, alpha=0.5)

    # Phase boundary lines
    ep1 = g0
    ep2 = g0 + gen_interv[1]
    ep3 = g0 + gen_interv[1] + gen_interv[2]
    for ep in [ep1, ep2, ep3]:
        ax[0].plot([ep, ep], [-1, 18.5], 'k-.', alpha=0.5)
        ax[1].plot([ep, ep], [40, 85.5], 'k-.', alpha=0.5)

    ax[0].set_ylim([-0.2, 4.0])
    ax[1].set_ylim([43, 92])
    ax[0].legend(loc='upper left', fontsize=16)
    ax[1].legend(loc='lower left', fontsize=16)
    ax[0].set_title('MSN\n\n$\\kappa$=' + my_kappa, fontsize=18)
    ax[1].set_title('GPi', fontsize=18)
    for a in ax:
        a.spines['right'].set_visible(False)
        a.spines['top'].set_visible(False)
        a.set_ylabel('Mean firing rate [Hz]', fontsize=16)
        a.tick_params(axis='both', which='major', labelsize=14)
    ax[1].set_xlabel('Episodes', fontsize=16)
    plt.subplots_adjust(hspace=0.35)
    if save:
        plt.savefig(fname + '_' + str(time.time()) + '.png', format='png', dpi=300)
    plt.show()


def _plot_stdp_2x2(burst_params, dip_params, suptitle, save, fname):
    """
    Helper for Figure 8 C/E/F: 2×2 STDP grid.
    burst_params/dip_params: (D1_Aplus, D1_Aminus, D2_Aplus, D2_Aminus)
    """
    dc  = 0.00056
    tau = 20.0
    dt  = np.arange(0., 40., 0.01)

    conditions = [
        ('D1-MSN (DA burst)', burst_params[0], burst_params[1], 'steelblue', 0, 0),
        ('D1-MSN (DA dip)',   dip_params[0],   dip_params[1],   'steelblue', 0, 1),
        ('D2-MSN (DA burst)', burst_params[2],  burst_params[3],  'crimson',   1, 0),
        ('D2-MSN (DA dip)',   dip_params[2],    dip_params[3],    'crimson',   1, 1),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(8, 8), sharex=True)
    if suptitle:
        fig.suptitle(suptitle, fontsize=16)

    for title, Aplus, Aminus, color, row, col in conditions:
        ax = axes[row][col]
        ax.plot(dt,      Aplus  * np.exp(-dt / tau), color=color, linewidth=3)
        ax.plot(-1.*dt, -Aminus * np.exp(-dt / tau), color=color, linewidth=3)
        ax.plot([0., 0.], [-dc, dc], 'k:')
        ax.plot([-40, 40], [0., 0.], 'k:')
        ax.set_title(title, fontsize=18)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks([0])
        ax.tick_params(axis='both', which='major', labelsize=15)
        if col == 0:
            ax.text(-42, dc * 0.6,  'A+', fontsize=16, va='center')
            ax.text(-42, -dc * 0.6, 'A-', fontsize=16, va='center')
        if row == 1:
            ax.set_xlabel('$\\Delta t$ (ms)', fontsize=17)

    plt.tight_layout()
    if save:
        plt.savefig(fname + '_' + str(time.time()) + '.png', format='png', dpi=300)
    plt.show()


def Figure_8_A(save=True):
    """Figure 8A: Parkinson model, λ=0.1, κ=4.0, single run (log0)."""
    _plot_individual_trajectory('run4_k_4_l_01_m_1_p_no_conv', 228, 0,
                                [100, 30, 100], 8, 'Parkinson λ=0.1', save, 'fig_8A')


def Figure_8_B(save=True):
    """Figure 8B: Parkinson model, λ=0.2, κ=4.0, single run (log2 — complete 230-episode run)."""
    _plot_individual_trajectory('run4_k_4_l_02_m_1_p_no_conv', 228, 2,
                                [100, 30, 100], 8, 'Parkinson λ=0.2', save, 'fig_8B')


def Figure_8_D(save=True):
    """Figure 8D: Schizophrenia model, λ=0.2, κ=4.0, single run (log0)."""
    _plot_individual_trajectory('run4_k_4_l_02_m_1_s', 228, 0,
                                [100, 30, 100], 8, 'Schizophrenia λ=0.2', save, 'fig_8D')


def Figure_8_C(save=True):
    """Figure 8C: STDP functions — Parkinson (reversed DA burst for D1, D2)."""
    dc = 0.00056
    burst = (-dc,       dc/4.0,   dc/3.0,  -dc*3.0/4.0)
    dip   = (-dc/40.0,  dc/40.0,  dc,      -dc)
    _plot_stdp_2x2(burst, dip, 'STDP functions (Parkinson)', save, 'fig_8C')


def Figure_8_E(save=True):
    """Figure 8E: STDP functions — Schizophrenia (reduced DA dip for D2)."""
    dc = 0.00056
    burst = (dc,        dc/4.0,   dc/3.0,   dc*3.0/4.0)
    dip   = (-dc/40.0,  dc/40.0,  dc/4.0,  -dc/4.0)
    _plot_stdp_2x2(burst, dip, 'STDP functions (Schizophrenia)', save, 'fig_8E')


def Figure_8_F(save=True):
    """Figure 8F: STDP functions — Normal conditions (same as Figure_3_D)."""
    dc = 0.00056
    burst = (dc,        dc/4.0,   dc/3.0,   dc*3.0/4.0)
    dip   = (-dc/40.0,  dc/40.0,  dc,       -dc)
    _plot_stdp_2x2(burst, dip, 'STDP functions (Normal Conditions)', save, 'fig_8F')


if __name__ == '__main__':
    Figure_4_A(save=False)
    Figure_1_E(save=False)
    Figure_1_F(save=False)
    Figure_1_G(save=False)
    Figure_3_C(save=False)
    Figure_3_D(save=False)
    Figure_4_D(save=False)
    Figure_5_A(save=False)
    Figure_6_AB(save=False)
    Figure_6_C(save=False)
    Figure_7_D(save=False)
    Figure_8_A(save=False)
    Figure_8_B(save=False)
    Figure_8_C(save=False)
    Figure_8_D(save=False)
    Figure_8_E(save=False)
    Figure_8_F(save=False)
