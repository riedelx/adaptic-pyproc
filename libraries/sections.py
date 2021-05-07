import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import utils
import math

def plotFunc(self,x_coordinates,y_coordinates,reverse,title='section layout', analysis_2d = True):
    fig = plt.figure(figsize = (6,4))
    ax = fig.add_subplot(111)
    ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
    if reverse:y_coordinates=[-i for i in y_coordinates]
    ax.plot(x_coordinates, y_coordinates,'k-')
    try:
        if self.analysis_2d:
            for i in self.reinf_sect:
                yTemp=i[2]
                bTemp=self.width(yTemp)
                bSpacing=bTemp/(i[0]+1)
                if reverse:yTemp=-yTemp
                for j in range(i[0]):
                    xTemp=-bTemp/2+j*bSpacing+bSpacing
                    ax.add_patch(plt.Circle((xTemp,yTemp), radius=i[1]/2, color='b', fill=False))
        else:
            for i in self.reinf_sect:
                yTemp = i[1]
                if reverse: yTemp=-yTemp
                for xTemp in i[2]:
                    ax.add_patch(plt.Circle((xTemp,yTemp), radius=i[0]/2, color='b', fill=False))
                    if xTemp != 0: ax.add_patch(plt.Circle((-xTemp,yTemp), radius=i[0]/2, color='b', fill=False))
        ax.add_artist(circ)
    except: pass
    if reverse:ax.plot(0,-self.centr,'r+',markersize=10,linewidth=8)
    else:ax.plot(0,self.centr,'r+',markersize=10,linewidth=8)
    ax.set_title(title)
    ax.set_aspect('equal', 'box')
    plt.show()

class rss:
    def __init__(self, ID, mat, b, d):
        self.ID = ID
        self.mat = mat
        self.b = b
        self.d = d
        self.area= b * d
        self.centr_top = d/2
        self.centr_bot = -d/2
        self.Iy = b*d**3/12
        self.Ix = b**3*d/12

    def plotting(self,title='section layout'):
            y_coordinates = [0,d,d,0,0]
            x_coordinates = [b/2,b/2,-b/2,-b/2,b/2]
            plotFunc(self,x_coordinates,y_coordinates,title=title)
    def width(self,x):
        if (x >= 0 and x <= self.d):
            b=self.b
        else:
            b=0
        return b

    def adaptic_print(self):
        # convert dictionary values of different dimensionality to a flatten list
        line = self.__dict__
        line = np.array(list(line.values())[0:4])
        return utils.str_joint(line)

class rccs:
    def __init__(self, ID, reinf_mat, unconf_mat, conf_mat, hc1, bc1, cover, links, reinf):
        self.ID = ID
        self.reinf_mat = reinf_mat
        self.unconf_mat = unconf_mat
        self.conf_mat = conf_mat
        self.hc1 = hc1
        self.hc2 = hc1 -2 * cover - links
        self.bc1 = bc1
        self.bc2 = bc1 -2 * cover - links
        # each layer of reinforcement given as: [no of bars, diameter, distance from the bottom fibre]
        # reinforcement given as a list [[layer_1], [layer_2],..., [layer_n]]
        # layers order starting from the furtherst layer from the bottom fibre
        reinf = pd.DataFrame(reinf).sort_values(by=[2],ascending=False).astype(int).values.tolist()
        self.reinf = [[int(math.pi*i[1]**2/4 * i[0]), i[2]] for i in reinf]
        self.reinf_sect = reinf # [no of bars, diameter, distance from the bottom fibre]
        self.area = hc1 * bc1
        self.cover = cover
        self.links = links # links diameter

    def adaptic_print(self):
        # convert dictionary values of different dimensionality to a flatten list
        dict0 = self.__dict__
        dict1 = np.array(list(dict0.values())[0:8])
        dict2 = np.array(list(dict0.values())[8])
        dict2 = [val for sublist in dict2 for val in sublist]
        line = np.hstack((dict1,dict2))
        return utils.str_joint(line)

class rcts:
    def __init__(self, ID, reinf_mat, unconf_mat, conf_mat, Df, Dw, Bf, Bw, cover, links, reinf_sect, analysis_2d = True):
        self.ID = ID
        self.reinf_mat = reinf_mat
        self.unconf_mat = unconf_mat
        self.conf_mat = conf_mat
        self.Df = Df
        self.Dw = Dw
        self.df = Df - cover * 2 - links
        self.db = Dw # + cover
        self.Bf = Bf
        self.Bw = Bw
        self.bf = Bf - 2 * cover - links
        self.bw = Bw - 2 * cover - links

        # reinforcement given as a list [[layer_1], [layer_2],..., [layer_n]]
        # for 2D each layer of reinforcement given as: [no of bars, diameter, distance from the bottom fibre]
        # for 2D each layer of reinforcement given as: [diameter, distance from the bottom fibre, [z1, z2, z3]]
        # after sorting layers order starting from the furtherst layer from the bottom fibre

        if analysis_2d:
            reinf_sect.sort(key = lambda x: x[2], reverse = True)
            self.reinf = [[int(math.pi*i[1]**2/4 * i[0]), i[2]] for i in reinf_sect] # [no of bars, diameter, distance from the bottom fibre]
        else:
            reinf_sect.sort(key = lambda x: x[1], reverse = True)
            self.reinf = []
            for i in reinf_sect:
                for j in i[2]:
                    self.reinf.append([int(math.pi*i[0]**2/4), i[1], j])
        self.reinf_sect = reinf_sect
        self.analysis_2d = analysis_2d
        self.cover = cover
        self.links = links # links diameter
        area_w = Dw * Bw
        area_f = Df * Bf
        self.area = area_w + area_f
        self.centr_top = int((Df*area_f/2+area_w*(Df+Dw/2))/self.area)
        self.centr_bot = -int(Df + Dw- self.centr_top)

    def plotting(self,title='section layout',reverse=False):
        self.h=self.Dw+self.Df
        self.centr = self.centr_top #int(self.h-(Df*area_f/2+area_w*(Df+Dw/2))/self.area)
        y_coordinates = [0,self.Dw+0,self.Dw+0,self.h,self.h,self.Dw+0,self.Dw+0,0,0]
        x_coordinates = [self.Bw/2,self.Bw/2,self.Bf/2,self.Bf/2,-self.Bf/2,-self.Bf/2,-self.Bw/2,-self.Bw/2,self.Bw/2]
        plotFunc(self,x_coordinates,y_coordinates,title=title,reverse=reverse, analysis_2d = self.analysis_2d)
    def width(self,x):
        if (x >= 0 and x <= self.Dw):
            b=self.Bw
        elif x<=self.Dw+self.Df:
            b=self.Bf
        else:
            b=0
        return b

    def adaptic_print(self):
        # convert dictionary values of different dimensionality to a flatten list
        dict0 = self.__dict__
        dict1 = np.array(list(dict0.values())[0:12])
        dict2 = np.array(list(dict0.values())[12])
        dict2 = [val for sublist in dict2 for val in sublist]
        line = np.hstack((dict1,dict2))
        return utils.str_joint(line)

class isec:
    def __init__(self, ID, mat, bf1, tf1, bf2, tf2, dw, tw): # ID, Df, Dw, Bf, Bw, ):
        self.ID = ID
        self.mat = mat
        self.bf1 = bf1
        self.tf1 = tf1
        self.bf2 = bf2
        self.tf2 = tf2
        self.dw = dw
        self.tw = tw
        area_w = dw * tw
        area_f1 = bf1 * tf1
        area_f2 = bf2 * tf2
        self.area = area_w + area_f1 + area_f2
        self.centr_top = int(((bf1*tf1)*tf1/2+(bf2*tf2)*(tf1+dw+tf2/2)+(dw*tw)*(tf1+dw/2))/self.area)
        self.centr_bot = -int(tf1+tf2+dw-self.centr_top)

    def plotting(self,title='section layout',reverse=False):
        self.h=self.tf1+self.tf2+self.dw
        self.centr = self.centr_top
        y_coordinates = [0,self.tf1,self.tf1,self.tf1+self.dw,self.tf1+self.dw,self.h,self.h,self.tf1+self.dw,self.tf1+self.dw,self.tf1,self.tf1,0,0]
        x_coordinates = [self.bf1/2,self.bf1/2,self.tw/2,self.tw/2,self.bf2/2,self.bf2/2,-self.bf2/2,-self.bf2/2,-self.tw/2,-self.tw/2,-self.bf1/2,-self.bf1/2,self.bf1/2]
        plotFunc(self,x_coordinates,y_coordinates,title=title,reverse=reverse)
    def width(self,x):
        if (x >= 0 and x <= self.tf1):
            b=self.bf1
        if (x >= 0 and x <= self.tf1+self.dw):
            b=self.tw
        elif x<=self.tf1+self.tf2+self.dw:
            b=self.bf2
        else:
            b=0
        return b

    def adaptic_print(self):
        # convert dictionary values of different dimensionality to a flatten list
        dict0 = self.__dict__
        dict1 = np.array(list(dict0.values())[0:8])
        return utils.str_joint(dict1)
