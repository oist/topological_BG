
################################
### creator: Kenji Doya ########
### 2019-11-19 #################
################################

import numpy as np
import math

# Point mass model
class pmm2:
    """Class for a 2D point mass model"""

    def __init__(self, mass=0.1, friction=0.2):
        """Create a new point mass
        default: mass=0.1kg, friction=0.2N/(m/s)"""
        self.m = mass
        self.mu = friction
        self.x = np.zeros(2)  # position
        self.v = np.zeros(2)  # velocity

    def move(self, u, dt):
        """Move by input u=[ux, uy] for dt"""
        if math.isnan(u[0]) or math.isnan(u[1]):
            print('nan value will be skipped ..')
        else:
            self.x += self.v*dt
            self.v += (np.array(u) - self.mu*self.v)/self.m*dt

# BG reinforcement learning agent
class QLa:
    """Class for a Q-learning agent"""
    
    def __init__(self, nact=4, rate=0.2, itemp=2, tdisc=0):
        """Create a new RL agent
        default: nact=4, rate=0.1, itemp=1, tdisc=0"""
        self.na = nact
        self.alpha = rate
        self.beta = itemp
        self.gamma = tdisc
        self.Q = np.zeros(self.na)  # action values
        self.action = 0  # chosen action
        self.delta = 0   # reward prediction error

    def select(self, opts=1):
        """Select an action from current options"""
        q = self.Q*opts  # mask by option vector       ### replace self.Q by 3 mean firing rates [MSN_L,MSN_C,MSN_R]
        p = np.exp(self.beta*q)
        p = p/sum(p)   #normalize
        # sample by multinoulli (categorical) distribution
        #print 'p:',p
        s = np.random.multinomial(1, p)
        #print 's: ',s
        self.action = list(s).index(1)  # find the index of 1   ###done by BG, GPI, find the index lowest rate channel index
        
        return self.action

    def update(self, r):
        """Update action value for reward r"""
        self.delta = r - self.Q[self.action]             ### release dopamine r.
        self.Q[self.action] += self.alpha*self.delta     ### record the new self.Q by 3 mean firing rates [MSN_L,MSN_C,MSN_R]
        #print self.Q

def reward(a):
    """reward function"""
    return (1.5 if a==1 else 0)  # 1 is the target  ### given by M2.
    #"""Reward for hand position x"""
    #xoff = x-xtarget[target]   # offset
    #return (1 if sum(xoff**2)<0.02**2 else 0)  # within 2cm


