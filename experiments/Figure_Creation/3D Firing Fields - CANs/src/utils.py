import numpy as np
import torch
from torch import nn
import torch.nn.functional as F
from src.twisted_attr import TwistedAttractor

class MultiAttractor:
    def __init__(self, m=3, n=20, dt=0.4):
        self.attr = []
        self.m = m
        for i in range(self.m):
            attr = TwistedAttractor(
                        n=n, 
                        sigma=2., 
                        dt=dt, 
                        v=0.1, 
                        cell_pos=None)
            
            self.attr.append(attr)

        self.reset()
    
    def reset(self, y=None):
        if y is None: self.y = torch.rand(self.m, 1, 400)*0.01
        else: self.y = y.view(self.m, 1, 400).clone()
    
    def __call__(self, v, b, dt=None):
        for i in range(self.m):
            self.y[i] = self.attr[i](self.y[i], v[i], b[i], dt)
            
        return self.y.view(1, self.m*400)



def get_path(T=1000, w=35, vmax=0.4):
    x = torch.zeros(T, 3)
    x[0] = w/2.

    v    = torch.zeros(T, 3)
    v[0] = torch.randn(3)
    v    = torch.clamp(v, -vmax, vmax)


    for t in range(1,T):
        a = torch.randn(3)*0.05
        for i in range(3):
            if x[t-1,i] <= 0:  a[i] =  0.2
            if x[t-1,i] >= w:  a[i] = -0.2

        v[t] = v[t-1] + a
        v[t] = torch.clamp(v[t], -vmax, vmax)
        x[t] = x[t-1] + v[t]
        
    return x, v