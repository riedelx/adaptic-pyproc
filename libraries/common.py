#============================================
# Common Functions
#============================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import log10, floor
from decimal import Decimal

def typename(x): return type(x).__name__
def convertFloat(number):
    if isinstance(number, float) or isinstance(number, int):
        if 1E5 < number or number < 1E-5:
            return '%.2E' % Decimal(number)
        elif isinstance(number, float):
            sig=3
            return str(round(number, sig-int(floor(log10(abs(number))))-1))
        else:
            return str(number)
    else:
        return str(number)
def str_joint(lst): return " ".join([convertFloat(i) for i in lst])
def data_frame(obj): return pd.DataFrame(obj.__dict__,index=[0])
def data_frame_alt(obj): return pd.DataFrame.from_dict(obj.__dict__)
def dict2lst(obj): return list(obj.__dict__.values())
def csv2np(fname):
    with open(fname) as f:
        # skip comment line manually using generator expression
        lines = (line for line in f if not line.startswith('#')) 
        FH = np.loadtxt(lines, delimiter='	', skiprows=1)
    return FH
 
def custom_plot():
    fig, ax = plt.subplots() 
    ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
    
def closestPoint(point, listTo, listFrom):
    index = next(x[0] for x in enumerate(listFrom) if x[1] > point)
    return np.array([float(listTo[index]), float(listFrom[index])])

def df_index(df,val,col_ID): return df.index[df[col_ID] == val].tolist()[0]

def df_value(df,val,col_ID0, col_ID1): 
    return df.loc[df_index(df,val,col_ID0)][col_ID1]

def name_gen(string,lst, start_ID = 1):
    return [string + x for x in np.arange(start_ID,start_ID + len(lst),1).astype(str)]
    
def create_ASTR(positive, negative = []):
    if negative == []:
        negative = np.negative(positive)
    ASTR = ['astr']
    for i in [positive, negative]:
        x1 = i[0][0]
        y1 = i[0][1]
        k1 = (y1)/(x1)
        x2 = i[1][0]
        y2 = i[1][1]
        k2 = (y2 - y1)/(x2 - x1)
        x3 = i[2][0]
        y3 = i[2][1]
        k3 = (y3 - y2)/(x3 - x2)
        ASTR.extend((k1,x1,k2,x2,k3))
    return ASTR
    
def ASTR_plot(lst, ult_p, ult_n, title, x_label, y_label, scaleX = 1, scaleY = 1):
    lst = lst[1:]
    positive = lst[0:5]
    negative = lst[5:10]
    positive_bool = True
    for i in [positive, negative]:
        x1 = i[1]
        k1 = i[0]
        y1 = k1 * x1
        x2 = i[3]
        k2 = i[2]
        y2 = y1 + k2 * (x2 - x1)
        k3 = i[4]
        if positive_bool:
            x3 = ult_p
            y3 = y2 + k3 * (x3 - x2)
            ASTR_pos = np.array([[0,0],[x1,y1],[x2,y2],[x3,y3]])
        else:
            x3 = ult_n
            y3 = y2 + k3 * (x3 - x2)
            ASTR_neg = np.flip(np.array([[0,0],[x1,y1],[x2,y2],[x3,y3]]), axis = 0)
        positive_bool = False
    ASTR_xy = np.vstack((ASTR_neg,ASTR_pos))
    fig, ax = plt.subplots() 
    ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
    plt.plot(ASTR_xy[:,0] * scaleX,ASTR_xy[:,1] * scaleY)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.show()

def setattrs(_self, **kwargs):
    for k,v in kwargs.items():
        setattr(_self, k, v)

class param2adap:
    def __init__(self, name, params):
        self.name = name
        self.params = params
        crv_type = []
        crv_param = []
        for i in params:
            if isinstance(i, str):
                crv_type.append(i)
            else:
                crv_type.append(i[0])
                crv_param.append([x for i,x in enumerate(i) if i!=0])
        self.crv_type = crv_type
        self.crv_param = [val for sublist in crv_param for val in sublist]
        self.printout = str_joint([name] + self.crv_type + self.crv_param)