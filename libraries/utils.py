#============================================
# Common Functions
#============================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import log10, floor
from decimal import Decimal

def checkSubset(list1, list2): #  Check if a nested list is a subset of another nested list
    if list1 == []: exist = False
    else:
        l1, l2 = list1[0], list2[0]
        exist = True
        for i in list2:
            if i not in list1:
                exist = False
    return exist

def removeRows(df,labels,colID='ID'):
    idxs = []
    for i in labels:
        idxs.append(df_index(df,i,colID))
    return df.drop(labels=idxs,axis=0).reset_index(drop=True)
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

def find_nearest(array, value): #numpy
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def closestPoint(point, listTo, listFrom):
    index = next(x[0] for x in enumerate(listFrom) if x[1] > point)
    return np.array([float(listTo[index]), float(listFrom[index])])

def df_index(df,val,col_ID): return df.index[df[col_ID] == val].tolist()[0]

def df_value(df,val,col_ID0, col_ID1):
    return df.loc[df_index(df,val,col_ID0)][col_ID1]

def name_gen(string,lst, start_ID = 1):
    return [string + x for x in np.arange(start_ID,start_ID + len(lst),1).astype(str)]

# disp1 = 1 # mm
# disp2 = 2 # mm
# disp3 = 3 # mm
# ultDispPos = disp3
# ultDispNeg = -disp3
# f1 = 5 # N
# f2 = 6 # N
# f3 = 7 # N
# negative=np.array([[-1,-0],[-5,-0],[-5,-0]])
# astr1 = utils.create_ASTR(np.array([[disp1, f1],[disp2, f2],[disp3, f3]]),negative=negative)
# utils.ASTR_plot(astr1, ultDispPos, ultDispNeg, 'ASTR curve','displacement [mm]', 'force [kN]', scaleX = 1, scaleY = 1E-3)

def create_ASTR(positive, negative = []):
    # create_ASTR(np.array([[disp1, f1],[disp2, f2],[disp3, f3]]),negative=np.array([[-disp1, -f1],[-disp2, -f2],[-disp3, -f3]]))
    # utils.create_ASTR(np.array([[disp1, f1],[disp2, f2],[disp3, f3]]),negative=np.array([['soft']]))
    if len(negative) == 0:
        negative = np.negative(positive)
    ASTR = ['astr']
    for i in [positive, negative]:
        if len(i) == 1 and i == ['soft']:
            k1 = 1
            x1 = -0
            k2 = 0
            x2 = -1
            k3 = 0
        else:
            x1 = i[0][0]
            y1 = i[0][1]
            k1 = (y1)/(x1)
            x2 = i[1][0]
            y2 = i[1][1]
            if x2==x1:
                k2=k1
            else:
                k2 = (y2 - y1)/(x2 - x1)
            x3 = i[2][0]
            y3 = i[2][1]
            if x2==x3:
                k3=k2
            else:
                k3 = (y3 - y2)/(x3 - x2)
        ASTR.extend((k1,x1,k2,x2,k3))
    return ASTR

def ASTR_plot(lst, ult_p, ult_n, title='astr', x_label='', y_label='', scaleX = 1, scaleY = 1,xlim=(None,None),plotSide='both',returnData = False):
    #lst = lst[1:]
    positive = lst[1:6]
    negative = lst[6:11]
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
            if ult_p >= x2:
                x3 = ult_p
            else:
                x3 = x2
            y3 = y2 + k3 * (x3 - x2)
            # check if the line crosses 0
            if y3/y2 < 0:
                y3, y4 = 0, 0
                x4 = x2-y2/k3
            else:
                x4, y4 = x3, y3
            ASTR_pos = np.array([[0,0],[x1,y1],[x2,y2],[x4,y4],[x3,y3]])
        else:
            if np.abs(ult_n) >= np.abs(x2):
                x3 = ult_n
            else:
                x3 = x2
            y3 = y2 + k3 * (x3 - x2)
            # check if the line crosses 0
            if y3/y2 < 0:
                y3 = 0
                y4 = 0
                x4 = x2-y2/k3
            else:
                x4, y4 = x3, y3
            ASTR_neg = np.flip(np.array([[0,0],[x1,y1],[x2,y2],[x4,y4],[x3,y3]]), axis = 0)
        positive_bool = False
    if plotSide=='both':
        ASTR_xy = np.vstack((ASTR_neg,ASTR_pos))
    elif plotSide=='positive':
        ASTR_xy = ASTR_pos
    elif plotSide=='negative':
        ASTR_xy = ASTR_neg
    if returnData: return ASTR_xy
    else:
        fig, ax = plt.subplots()
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
        plt.plot(ASTR_xy[:,0] * scaleX,ASTR_xy[:,1] * scaleY)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xlim(xlim)
        plt.title(title)
        plt.show()

def astr_gap(gapP,gapN,S1=2.00E+09,S2=1.00E+09,plotting=True,xlim=[20,-20],title='astr gap curve',xlabel='displacement [mm]',ylabel='force [kN]', scaleX = 1, scaleY = 1):
    # S1, S2 are the stiffnesses after closing the gap, S1 > S2
    # gapP > 0, gapN < 0 or rigid
    # utils.astr_gap(9999,-0.001,S1=2E+9,S2=1E+9,plotting=False)
    k,gap=[],[]
    for i in [gapP,gapN]:
        if i=='rigid':
            k.append(1.00E+09)
            gap.append(0.1)
        else:
            k.append(0.1)
            gap.append(np.abs(i))
    curve=['astr']+[S1, 0.00E+00, k[0], gap[0], S2]+[S1, 0.00E+00, k[1], -gap[1], S2]
    if plotting:
             ASTR_plot(curve, xlim[0], xlim[1], title,xlabel, ylabel, scaleX = scaleX, scaleY = scaleY)
    return curve

def setattrs(_self, **kwargs):
    for k,v in kwargs.items():
        setattr(_self, k, v)

def centroidX(points,discr=100):
    #points=[[x1,y1],[x2,y2]]
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    points=np.array([x,y])
    area=np.trapz([y[0],y[-1]], x=[x[0],x[-1]])
    L=x[-1]-x[0]
    L_incr=L/discr
    moment=0
    for i in range(discr):
        if i == discr-1:
            x2=x[-1]
        else:
            x2=(i+1)*L_incr
        y2=findExactPoint(points, x2,limY=False)[1]
        x1=(i)*L_incr
        y1=findExactPoint(points, x1,limY=False)[1]
        area_temp=np.trapz([y1,y2], x=[x1,x2])
        moment+=area_temp*(x2-L_incr/2)
    return moment/area

# #=====================================================#
# Intersection of two curves
"""
Sukhbinder
5 April 2017
Based on:
"""
def _rect_inter_inner(x1,x2):
    n1=x1.shape[0]-1
    n2=x2.shape[0]-1
    X1=np.c_[x1[:-1],x1[1:]]
    X2=np.c_[x2[:-1],x2[1:]]
    S1=np.tile(X1.min(axis=1),(n2,1)).T
    S2=np.tile(X2.max(axis=1),(n1,1))
    S3=np.tile(X1.max(axis=1),(n2,1)).T
    S4=np.tile(X2.min(axis=1),(n1,1))
    return S1,S2,S3,S4

def _rectangle_intersection_(x1,y1,x2,y2):
    S1,S2,S3,S4=_rect_inter_inner(x1,x2)
    S5,S6,S7,S8=_rect_inter_inner(y1,y2)

    C1=np.less_equal(S1,S2)
    C2=np.greater_equal(S3,S4)
    C3=np.less_equal(S5,S6)
    C4=np.greater_equal(S7,S8)

    ii,jj=np.nonzero(C1 & C2 & C3 & C4)
    return ii,jj

def name_reset(df1,df2,prefix="nst"):
    df1=df1.drop(["name"], axis=1)
    df=pd.concat([df1,df2])
    df=df.drop_duplicates(inplace=False,subset='ID')
    dfName = utils.name_gen(prefix,df)
    df.insert(0, "name", dfName, True)
    return df.reset_index(drop=True)

def intersection(x1,y1,x2,y2):
    """
INTERSECTIONS Intersections of curves.
   Computes the (x,y) locations where two curves intersect.  The curves
   can be broken with NaNs or have vertical segments.
usage:
x,y=intersection(x1,y1,x2,y2)
    Example:
    a, b = 1, 2
    phi = np.linspace(3, 10, 100)
    x1 = a*phi - b*np.sin(phi)
    y1 = a - b*np.cos(phi)
    x2=phi
    y2=np.sin(phi)+2
    x,y=intersection(x1,y1,x2,y2)
    plt.plot(x1,y1,c='r')
    plt.plot(x2,y2,c='g')
    plt.plot(x,y,'*k')
    plt.show()
    """
    x1,y1,x2,y2 = np.array(x1),np.array(y1),np.array(x2),np.array(y2)
    ii,jj=_rectangle_intersection_(x1,y1,x2,y2)
    n=len(ii)

    dxy1=np.diff(np.c_[x1,y1],axis=0)
    dxy2=np.diff(np.c_[x2,y2],axis=0)

    T=np.zeros((4,n))
    AA=np.zeros((4,4,n))
    AA[0:2,2,:]=-1
    AA[2:4,3,:]=-1
    AA[0::2,0,:]=dxy1[ii,:].T
    AA[1::2,1,:]=dxy2[jj,:].T

    BB=np.zeros((4,n))
    BB[0,:]=-x1[ii].ravel()
    BB[1,:]=-x2[jj].ravel()
    BB[2,:]=-y1[ii].ravel()
    BB[3,:]=-y2[jj].ravel()

    for i in range(n):
        try:
            T[:,i]=np.linalg.solve(AA[:,:,i],BB[:,i])
        except:
            T[:,i]=np.NaN


    in_range= (T[0,:] >=0) & (T[1,:] >=0) & (T[0,:] <=1) & (T[1,:] <=1)

    xy0=T[2:,in_range]
    xy0=xy0.T
    return xy0[:,0],xy0[:,1]

def findExactPoint(curve1, coordinate,limY=True, multiple=False):
    if limY:
        curve2=np.array([[-9E99,np.inf],[coordinate,coordinate]])
    else:
        curve2=np.array([[coordinate,coordinate],[-9E99,np.inf]])
    tupl = (intersection(curve1[0], curve1[1], curve2[0], curve2[1]))
    try:
        return [float(tupl[0]),float(tupl[1])]
    except:
        if multiple:
            return (tupl[0]),(tupl[1])
        else:
            return [float(tupl[0][0]),float(tupl[1][0])]
# #=====================================================#

class param2adap:
    def __init__(self, name, params,ID=''):
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
        if ID=='': ID=''
        else: ID=' # '+ID
        self.printout = str_joint([name] + self.crv_type + self.crv_param)+ID

# Maths and geometry
import math
def dot(v,w):
    x,y,z = v
    X,Y,Z = w
    return x*X + y*Y + z*Z

def length(v):
    x,y,z = v
    return math.sqrt(x*x + y*y + z*z)

def vector(b,e):
    x,y,z = b
    X,Y,Z = e
    return (X-x, Y-y, Z-z)

def unit(v):
    x,y,z = v
    mag = length(v)
    return (x/mag, y/mag, z/mag)

def distance(p0,p1):
    return length(vector(p0,p1))

def scale(v,sc):
    x,y,z = v
    return (x * sc, y * sc, z * sc)

def add(v,w):
    x,y,z = v
    X,Y,Z = w
    return (x+X, y+Y, z+Z)


# Given a line with coordinates 'start' and 'end' and the
# coordinates of a point 'pnt' the proc returns the shortest
# distance from pnt to the line and the coordinates of the
# nearest point on the line.
#
# 1  Convert the line segment to a vector ('line_vec').
# 2  Create a vector connecting start to pnt ('pnt_vec').
# 3  Find the length of the line vector ('line_len').
# 4  Convert line_vec to a unit vector ('line_unitvec').
# 5  Scale pnt_vec by line_len ('pnt_vec_scaled').
# 6  Get the dot product of line_unitvec and pnt_vec_scaled ('t').
# 7  Ensure t is in the range 0 to 1.
# 8  Use t to get the nearest location on the line to the end
#    of vector pnt_vec_scaled ('nearest').
# 9  Calculate the distance from nearest to pnt_vec_scaled.
# 10 Translate nearest back to the start/end line.
# Malcolm Kesson 16 Dec 2012

def pnt2line(pnt, start, end):
    line_vec = vector(start, end)
    pnt_vec = vector(start, pnt)
    line_len = length(line_vec)
    line_unitvec = unit(line_vec)
    pnt_vec_scaled = scale(pnt_vec, 1.0/line_len)
    t = dot(line_unitvec, pnt_vec_scaled)
    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0
    nearest = scale(line_vec, t)
    dist = distance(nearest, pnt_vec)
    nearest = add(nearest, start)
    return (dist, nearest)

# Rather than using a for loop, you can vectorize these operations and get much better performance.
# Here is my solution that allows you to compute the distance from a single point to multiple line segments with vectorized computation.
def lineseg_dists(p, a, b):
    """Cartesian distance from point to line segment

    Edited to support arguments as series, from:
    https://stackoverflow.com/a/54442561/11208892

    Args:
        - p: np.array of single point, shape (2,) or 2D array, shape (x, 2)
        - a: np.array of shape (x, 2)
        - b: np.array of shape (x, 2)
    """
    # normalized tangent vectors
    d_ba = b - a
    d = np.divide(d_ba, (np.hypot(d_ba[:, 0], d_ba[:, 1])
                           .reshape(-1, 1)))

    # signed parallel distance components
    # rowwise dot products of 2D vectors
    s = np.multiply(a - p, d).sum(axis=1)
    t = np.multiply(p - b, d).sum(axis=1)

    # clamped parallel distance
    h = np.maximum.reduce([s, t, np.zeros(len(s))])

    # perpendicular distance component
    # rowwise cross products of 2D vectors
    d_pa = p - a
    c = d_pa[:, 0] * d[:, 1] - d_pa[:, 1] * d[:, 0]

    return np.hypot(h, c)
