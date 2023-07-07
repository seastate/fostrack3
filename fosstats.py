import matplotlib
matplotlib.use('TKAgg')
from matplotlib.patches import Ellipse
from matplotlib import pyplot as plt
plt.ion()
import numpy as np
from math import *


def readPathFile(pathfile):
    pathlist = []
    pf = open(pathfile,'r')
    while True:
        # read in the path number
        line = pf.readline()
        if not line:
            break
        pathnum = int(line)
        print(f'pathnum = {pathnum}')
        line = pf.readline()
        print(line)
        print( line.split(' ')[:-1])
        frames = [int(frmstr) for frmstr in line.split(' ')[:-1]]
        line = pf.readline()
        xs = [float(xstr) for xstr in line.split(' ')[:-1]]
        line = pf.readline()
        ys = [float(ystr) for ystr in line.split(' ')[:-1]]
        pathlist.append(SPath(pathnum=pathnum,frames=frames,xs=xs,ys=ys))
    pf.close()
    return pathlist


class SPath():
    def __init__(self,pathnum=None,frames=[],xs=[],ys=[]):
        self.pathnum = pathnum
        self.frames = frames
        self.xs = xs
        self.ys = ys
        self.pathlen = len(self.frames)

    def spans_frame(self,frame):
        """Return True iff this path spans frame"""
        if min(self.frames)<=frame and max(self.frames)>=frame:
            return True
        else:
            return False
        
    def stats(self):
        self.duration = self.frames[-1]-self.frames[0]
        self.net_disp = sqrt((self.xs[0]-self.xs[-1])**2+(self.ys[0]-self.ys[-1])**2)
        self.tot_disp = 0.
        for i in range(self.pathlen-1):
            self.tot_disp += sqrt((self.xs[i]-self.xs[i+1])**2+(self.ys[i]-self.ys[i+1])**2)
        
    def plot(self,ax):
        ax.plot(self.xs,self.ys)
        ax.annotate(str(self.pathnum),xy=(self.xs[0],self.ys[0]),xytext=(-3.,3.),textcoords='offset points')
        












