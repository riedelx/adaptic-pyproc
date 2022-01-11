import numpy as np
import pandas as pd
import utils
import matplotlib.pyplot as plt

# class con1:
#     def __init__(self, ID, fc1, length, fc2_factor = 0.05,epsilon_2t='',characteristic = True, ft_factor=1, Ec2 = '',Et2 = '',strain_prec=5,rev17=True,GfFactor=1,Qs=0.005,fy=500,tensionStiff=False):
#         self.resid_str = fc2_factor
#         self.ID = ID
#         self.fc1 = fc1
#         self.fc2 = self.resid_str * fc1
#         self.Qs=Qs # Qst - ratio of volume of transverse reinforcement to volume of concrete core
#         self.fy=fy # fy = yield strength of steel
#         self.length = length
#         self.density = 2400 / 10**9 # 2400 kg/m3
#         self.rev17=rev17
#         if characteristic:
#             self.fcm = fc1+8
#         else:
#             self.fcm = fc1
#         if fc1 <= 50:
#             self.ft = ft_factor*0.3 * self.fcm ** (2/3)
#         else:
#             self.ft = ft_factor*2.12*np.log(1+0.1*self.fcm)
#         if tensionStiff: self.ft = self.ft / 2
#         self.Gf = 73 * self.fcm**0.18/1000*GfFactor
#         self.Ec0 = int(21500*(self.fcm/10)**(1/3))
#         self.poisson = 0.2
#         self.Gc = 250 * self.Gf
#         self.epsilon_1c = 5 * fc1 / self.Ec0 /3
#         self.Ec1 = fc1 / self.epsilon_1c
#         self.alphaC = min(max(0,(self.Ec0 - self.Ec1)/self.Ec1),1)
#         if Ec2 != '':
#             self.Ec2 = Ec2
#             self.epsilon_2c = (self.fc1-self.fc2)/-Ec2+self.epsilon_1c
#         else:
#             # epsilon_2c - compressive fracture energy method
#             # self.Gc = round(250 * self.Gf, 1)
#             # # self.epsilon_2c = round((self.fc1-self.fc2)/-Ec2+self.epsilon_1c, strain_prec)
#             # self.epsilon_2c = round(self.epsilon_1c + 3 * self.Gc / (2 * length * fc1), strain_prec)
#
#             # epsilon_2c - Scott et al. (1982)
#             self.epsilon_2c = 0.004 + 0.9*self.Qs*self.fy/300
#             self.Ec2 = -(self.fc1-self.fc2)/(self.epsilon_2c - self.epsilon_1c) # secant compressive softening stiffness
#         self.Et1 = self.Ec0
#         self.epsilon_1t = self.ft / self.Et1
#         if epsilon_2t!='': self.epsilon_2t = epsilon_2t
#         elif tensionStiff:
#             self.epsilon_2t = 0.001
#         else:
#             alphat = -1
#             # epsilon_2t
#             area_f = self.Gf/self.length # fracture energy area
#             area_f_soft = area_f - self.epsilon_1t*self.ft/2 # area under the softening curve
#             # Method 1 -# Figure 1 page p 17 in RTD 2010
#             # eps_u = 2*area_f_soft/self.ft
#             # E0 = -self.ft/(eps_u-self.epsilon_1t) # tangent stiffness at epsilon_1t
#             # E1 = 0 # tangent stiffness at epsilon_2t
#             # self.epsilon_2t = max((self.epsilon_1t*E0+self.epsilon_1t*E1-2*self.ft)/(E0+E1),self.epsilon_1t)
#             # Method 2
#             # epsilon_2t using area under parabola curve
#             if area_f_soft > 0: self.epsilon_2t = (self.epsilon_1t*self.ft*alphat+3*self.epsilon_1t*self.ft+6*area_f_soft)/(self.ft*(alphat+3))
#             # self.epsilon_2t = max(self.epsilon_1t+3*area_f_soft/self.ft,self.epsilon_1t)
#             else: self.epsilon_2t = self.epsilon_1t
#         self.epsilon_2t = self.epsilon_2t
#         if Et2 != '':
#             self.Et2 = Et2
#             self.epsilon_2t = (self.ft)/-Et2+self.epsilon_1t
#         else:
#             self.Et2 = - self.ft /(self.epsilon_2t - self.epsilon_1t)
#
#     def adaptic_print(self,rawNr=False):
#         if rawNr: return self.ID,'con1', self.Ec1, self.fc1, self.Ec2, self.fc2, self.Et1,self.ft, self.Et2
#         elif self.rev17: line = utils.str_joint([self.ID,'con1', self.Ec1, self.fc1, self.Ec2, self.fc2, self.Et1,
#                                       self.ft, self.Et2, -self.alphaC, -1])
#         else: line = utils.str_joint([self.ID,'con1', self.Ec1, self.fc1, self.Ec2, self.fc2, self.Et1,
#                                       self.ft, self.Et2, self.alphaC])
#         return line

class con1:
    def __init__(self, ID, fc1, fcm, length, fc2_factor = 0.05, epsilon_2c=0.01, tensionStiff=False, customMaterial=False): # , Ec2_factor=1
        self.resid_str = fc2_factor
        self.ID = ID
        self.fc1 = fc1
        self.fc2 = self.resid_str * fc1
        self.length = length
        self.density = 2400 / 10**9 # 2400 kg/m3
        self.fcm = fcm
        if fc1 <= 50:
            self.ft = 0.3 * self.fcm ** (2/3)
        else:
            self.ft = 2.12*np.log(1+0.1*self.fcm)
        if tensionStiff: self.ft = self.ft / 2
        self.Gf = 73 * self.fcm**0.18/1000
        self.Ec0 = 21500*(self.fcm/10)**(1/3)
        # self.poisson = 0.2

        # COMPRESSIVE HARDENING
        # Method 1 - find alpha
        # self.epsilon_1c = 5 * fc1 / self.Ec0 /3
        # self.Ec1 = fc1 / self.epsilon_1c
        # self.alphaC = -min(max(0,(self.Ec0 - self.Ec1)/self.Ec1),1)

        # Method 2 - find Ec
        self.alphaC = -0.99
        if customMaterial=='plastic': self.alphaC = -0.99
        self.Ec1 = self.Ec0 / (np.abs(self.alphaC)+1)
        self.epsilon_1c =  fc1 / self.Ec1
        self.Ec1 = fc1 / self.epsilon_1c

        # COMPRESSIVE SOFTENING
        # Method 1 - compressive fracture energy method
        self.Gc = 250 * self.Gf * 4
        # # self.epsilon_2c = (self.fc1-self.fc2)/-self.Ec2+self.epsilon_1c
        # # self.epsilon_2c = self.epsilon_1c + 3 * self.Gc / (2 * length * fc1)
        self.epsilon_2c = (self.Gc/self.length-self.epsilon_1c*self.fc1)/((self.fc1+self.fc2)/2)+self.epsilon_1c

        # Method 2 - Scott et al. (1982)
        # self.epsilon_2c = 0.004 + 0.9*self.Qs*self.fy/300

        # Method 3
        # self.epsilon_2c = epsilon_2c

        # Finish softening
        self.Ec2 = -(self.fc1-self.fc2)/(self.epsilon_2c - self.epsilon_1c) # *Ec2_factor # secant compressive softening stiffness

        if customMaterial=='plastic':
            customMaterial=True
            self.Ec2 = -1.00E+01
            self.epsilon_2c = (self.fc1-self.fc2)/-self.Ec2+self.epsilon_1c

        # TENSION
        self.Et1 = self.Ec0
        self.epsilon_1t = self.ft / self.Et1
        self.alphaT = -1
        if tensionStiff:
            self.epsilon_2t = 0.001
        else:
            # epsilon_2t
            area_f = self.Gf/self.length # fracture energy area
            area_f_soft = area_f - self.epsilon_1t*self.ft/2 # area under the softening curve
            # Method 1 -# Figure 1 page p 17 in RTD 2010
            # eps_u = 2*area_f_soft/self.ft
            # E0 = -self.ft/(eps_u-self.epsilon_1t) # tangent stiffness at epsilon_1t
            # E1 = 0 # tangent stiffness at epsilon_2t
            # self.epsilon_2t = max((self.epsilon_1t*E0+self.epsilon_1t*E1-2*self.ft)/(E0+E1),self.epsilon_1t)
            # Method 2
            # epsilon_2t using area under parabola curve
            if area_f_soft > 0: self.epsilon_2t = (self.epsilon_1t*self.ft*self.alphaT+3*self.epsilon_1t*self.ft+6*area_f_soft)/(self.ft*(self.alphaT+3))
            # self.epsilon_2t = max(self.epsilon_1t+3*area_f_soft/self.ft,self.epsilon_1t)
            else: # area_f smaller than ft
                self.epsilon_1t = 2 * area_f / self.Et1
                self.epsilon_2t = 1.1*self.epsilon_1t
        self.epsilon_2t = self.epsilon_2t
        self.Et2 = - self.ft /(self.epsilon_2t - self.epsilon_1t)

        if customMaterial:
            self.Et2 = -1E+1
            self.ft = 0.1


        self.adaptic = [self.Ec1,self.fc1,self.Ec2,self.fc2,self.Et1,self.ft,self.Et2,self.alphaC,self.alphaT]

    def adaptic_print(self,rawNr=False):
        if not rawNr: return utils.str_joint([self.ID,'con1', self.Ec1, self.fc1, self.Ec2, self.fc2, self.Et1, self.ft, self.Et2, self.alphaC, self.alphaT])
        else: return [self.ID,'con1', self.Ec1, self.fc1, self.Ec2, self.fc2, self.Et1, self.ft, self.Et2, self.alphaC, self.alphaT]

    def data_frame(self):
        data = np.array([[self.ID, self.length, self.fc1, self.fc2, self.ft, self.Ec0, self.Ec1,
                          self.Ec2, self.Et1, self.Et2, self.Gf, self.Gc, self.epsilon_1c,
                          self.epsilon_2c, self.epsilon_1t, self.epsilon_2t,  self.alphaC]])
        df = pd.DataFrame(data,index=data[:,0])
        df.columns = ['ID', '$$h[mm]$$', '$$f_{c1}[MPa]$$','$$f_{c2}[MPa]$$', '$$f_{t}[MPa]$$',
                      '$$E_{c0}[MPa]$$','$$E_{c1}[MPa]$$','$$E_{c2}[MPa]$$','$$E_{t1}[MPa]$$',
                      '$$E_{t2}[MPa]$$','$$G_{f}[N/mm]$$','$$G_{c}[N/mm]$$','$$e_{c1}$$',
                      '$$e_{c2}$$','$$e_{t1}$$', '$$e_{t2}$$', '$$alpha_{c}$$']
        return df

class con2(con1):
    def __init__(self, ID, fc1, length, epsilon_t2 = 0.001, fc2_factor = 0.1, ft_factor = 1, characteristic = True, Ec2 = '',Et2 = '',strain_prec=5,k=1):
        super(con2, self).__init__(ID, fc1, length, epsilon_t2 = epsilon_t2, fc2_factor = fc2_factor, ft_factor = ft_factor, characteristic = characteristic, Ec2 = Ec2,Et2 = Et2,strain_prec=strain_prec)
        self.k=k
    def adaptic_print(self):
        line = utils.str_joint([self.ID,'con2', self.fc1, self.ft, self.epsilon_1c, self.k])
        return line


class stl1_el:
    def __init__(self, ID, E1, fy):
        self.ID = ID
        self.E1 = E1
        self.fy = fy
        self.density = 8050 / 10**9 # 8050 kg/m3
        self.epsilon_y = round(fy / E1,4)

    def adaptic_print(self):
        return utils.str_joint([self.ID,'stl1', self.E1, self.fy])

class stl1:
    def __init__(self, ID, E1, fy, fu, epsilon_u,softening=False,E3=5):
        self.ID = ID
        self.E1 = E1
        self.fy = fy
        self.fu = fu
        self.density = 8050 / 10**9 # 8050 kg/m3
        self.epsilon_u = round(epsilon_u,3)
        self.epsilon_y = round(fy / E1,4)
        self.E2 = round((fu - fy) / (epsilon_u - self.epsilon_y),1)
        self.E3=E3
        self.mu = round(self.E2 / E1,7)
        self.softening=softening

    def adaptic_print(self):
        if self.softening: line = utils.str_joint([self.ID,'stl1', self.E1, self.fy, self.mu, self.epsilon_u, self.E3])
        else: line = utils.str_joint([self.ID,'stl1', self.E1, self.fy, self.mu])
        return line

    def data_frame(self,alpha=''):
        data = np.array([[self.ID, self.E1, self.E2, self.fy, self.fu, self.epsilon_y, self.epsilon_u,
                          self.mu]])
        df = pd.DataFrame(data,index=data[:,0])
        df.columns = ['ID', '$$E_{1}[MPa]$$', '$$E_{2}[MPa]$$', '$$f_{y}[MPa]$$', '$$f_{u}[MPa]$$',
                      '$$e_{y}$$','$$e_{u}$$','$$mu$$']
        return df

    def stress_df(self):
        df = pd.DataFrame([[0,0],[self.epsilon_y,self.fy],[self.epsilon_u,self.fu]],columns=['strain','stress'])
        return df

    def stress(self,strain):
        if 0<=strain<=self.epsilon_y: return strain*self.E1
        elif self.epsilon_y<strain<=self.epsilon_u: return self.fy+(strain-self.epsilon_y)*self.E2
        elif 0>=strain>=-self.epsilon_y: return strain*self.E1
        elif -self.epsilon_y>strain>=-self.epsilon_u: return -self.fy+(strain-self.epsilon_y)*self.E2
        else: return 0

class bond:
    def __init__(self,ID,c,f_cm,ft,L,dia,n_bars,case=1,redFact=1):
        # case - refer to Table 6.1-1 MC2010:
        # 0 - Marti
        # 1 - Pull-out, good bond
        # 2 - Pull-out, all other bond cond
        # 3 - Splitting, good bond cond, unconfined
        # 4 - Splitting, good bond cond, stirrups
        # 5 - Splitting, all other bond cond, unconfined
        # 6 - Splitting, all other bond cond, stirrups
        # 7 - Debonded
        self.ID=ID
        self.case=case
        self.f_cm = f_cm # MPa
        self.c = c
        self.L=L
        self.dia=dia
        self.n_bars=n_bars
        self.redFact=redFact

    def slip2tau(self,s):
        if 0 <= s <= self.s_1:
            tau = self.tau_max*(s/self.s_1_original)**self.alpha
        elif self.s_1 <= s <= self.s_2:
            tau = self.tau_bu_split
        elif self.s_2 <= s <= self.s_3:
            tau = self.tau_bu_split - (self.tau_bu_split-self.tau_bf)*(s-self.s_2)/(self.s_3-self.s_2)
        else:
            tau = self.tau_bf
        return tau

    def create_ASTR(self):
        if self.case == 1 or self.case == 2 or self.case == 4 or self.case == 6:
            disp1,disp2,disp3,f1,f2,f3 = self.s_1,self.s_3,self.s_3*1.1,self.tau_bu_split,self.tau_bf,self.tau_bf
            self.astr_stress = utils.create_ASTR(np.array([[disp1, f1],[disp2, f2],[disp3, f3]]),negative=np.array([[-disp1, -f1],[-disp2, -f2],[-disp3, -f3]]))
        else:
            disp1,disp2,disp3,f1,f2,f3 = 0.2*self.s_1,self.s_1,self.s_3,self.slip2tau(0.2*self.s_1),self.tau_bu_split,self.tau_bf
            self.astr_stress = utils.create_ASTR(np.array([[disp1, f1],[disp2, f2],[disp3, f3]]),negative=np.array([[-disp1, -f1],[-disp2, -f2],[-disp3, -f3]]))

        # self.ft=ft # only for Marti
        # if self.case ==0: # Marti
        #     slip = 2 #dia*500**2/(8*210000*2*self.ft) # Alogla slip p 148
        #     self.tau_max = 2*self.ft * self.redFact
        #     self.s_1 = slip/2 # mm
        #     self.s_2 = slip # mm
        #     self.s_3 = slip # mm
        #     self.tau_bf = self.ft * self.redFact
        #     disp1,disp2,disp3,f1,f2,f3 = self.s_1,self.s_3,self.s_3+1,self.tau_max,self.tau_bf,self.tau_bf
        #     self.astr_stress = utils.create_ASTR(np.array([[disp1, f1],[disp2, f2],[disp3, f3]]),negative=np.array([[-disp1, -f1],[-disp2, -f2],[-disp3, -f3]]))
    def generateBond(self):
        if self.case ==1:
            self.tau_max = 2.5 * self.f_cm**0.5 * self.redFact
            self.tau_bu_split = self.tau_max
            self.s_1 = 1 # mm
            self.s_1_original = self.s_1
            self.s_2 = 2 # mm
            self.s_3 = self.c # mm, clear distance between ribs
            self.alpha = 0.4
            self.tau_bf = 0.4 * self.tau_max
            self.create_ASTR()
        elif self.case ==2:
            self.tau_max = 2.5 * self.f_cm**0.5 * self.redFact
            self.tau_bu_split = self.tau_max
            self.s_1 = 1.8 # mm
            self.s_1_original = self.s_1
            self.s_2 = 3.6 # mm
            self.s_3 = self.c # mm, clear distance between ribs
            self.alpha = 0.4
            self.tau_bf = 0.4 * self.tau_max
            self.create_ASTR()
        elif self.case ==3:
            self.tau_max = 2.5 * self.f_cm**0.5
            self.tau_bu_split = 7*(self.f_cm/25)**0.25 * self.redFact
            self.alpha = 0.4
            self.tau_bf = 0
            self.s_1_original = 1
            x = np.linspace(0,self.s_1_original,1000)
            y = [self.tau_max*(s/self.s_1_original)**self.alpha for s in x]
            self.s_1 = utils.findExactPoint([x,y], self.tau_bu_split, limY=True)[0]
            self.s_2 = self.s_1 # mm
            self.s_3 = 1.2 * self.s_1 # mm
            self.create_ASTR()
        elif self.case ==4:
            self.tau_max = 2.5 * self.f_cm**0.5
            self.tau_bu_split = 8*(self.f_cm/25)**0.25 * self.redFact
            self.alpha = 0.4
            self.tau_bf = 0.4 * self.tau_bu_split
            self.s_1_original = 1
            x = np.linspace(0,self.s_1_original,1000)
            y = [self.tau_max*(s/self.s_1_original)**self.alpha for s in x]
            self.s_1 = utils.findExactPoint([x,y], self.tau_bu_split, limY=True)[0]
            self.s_2 = self.s_1 # mm
            self.s_3 = 0.5 * self.c # mm
            self.create_ASTR()
        elif self.case ==5:
            self.tau_max = 1.25 * self.f_cm**0.5
            self.tau_bu_split = 5*(self.f_cm/25)**0.25 * self.redFact
            self.alpha = 0.4
            self.tau_bf = 0
            self.s_1_original = 1.8
            x = np.linspace(0,self.s_1_original,1000)
            y = [self.tau_max*(s/self.s_1_original)**self.alpha for s in x]
            self.s_1 = utils.findExactPoint([x,y], self.tau_bu_split, limY=True)[0]
            self.s_2 = self.s_1 # mm
            self.s_3 = 1.2 * self.s_1 # mm
            self.create_ASTR()
        elif self.case ==6:
            self.tau_max = 1.25 * self.f_cm**0.5
            self.tau_bu_split = 5.5*(self.f_cm/25)**0.25 * self.redFact
            self.alpha = 0.4
            self.tau_bf = 0.4 * self.tau_bu_split
            self.s_1_original = 1.8
            x = np.linspace(0,self.s_1_original,1000)
            y = [self.tau_max*(s/self.s_1_original)**self.alpha for s in x]
            self.s_1 = utils.findExactPoint([x,y], self.tau_bu_split, limY=True)[0]
            self.s_2 = self.s_1 # mm
            self.s_3 = 0.5 * self.c # mm
            self.create_ASTR()
        else: print('Case error: {}'.format(self.case))
    def force2tau(self,force):
        return list(np.array(force)/(self.n_bars*self.dia*np.pi*self.L))
    def tau2force(self,tau):
        return list(np.array(tau)*(self.n_bars*self.dia*np.pi*self.L))
    def curve(self,stop=9,num=1000,title='bond stress–slip relationship',astrPlot=True,label1='MC2010',label2='ASTR', returnData = False):
        fig = plt.figure(figsize = (6,4))
        ax = fig.add_subplot(111)
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
        x=np.linspace(0, stop, num)
        y=[self.slip2tau(i) for i in x]
        ax.plot(x,y,'-', linewidth=2, markersize=5,label=label1)
        if astrPlot:
            ASTR_xy = utils.ASTR_plot(self.astr_stress, stop, -stop, plotSide='positive', returnData = True)
            plt.plot(ASTR_xy[:,0],ASTR_xy[:,1])
            ax.legend()
        ax.set_title(title)
        ax.set_xlabel('Slip [mm]')
        ax.set_ylabel('Stress [MPa]')
        ax.set_xlim(0,None)
        ax.set_ylim(0,None)
        plt.show()
        if returnData: return [x,y]
    def astr_force(self,stop=0):
        if stop == 0: stop=self.s_3*1.1
        ASTR_x = utils.ASTR_plot(self.astr_stress, stop, 0, plotSide='positive', returnData = True)[:,0][1:]
        ASTR_y = utils.ASTR_plot(self.astr_stress, stop, 0, plotSide='positive', returnData = True)[:,1][1:]
        disp1,disp2,disp3,disp4 = ASTR_x
        f1,f2,f3,f4 = self.tau2force(ASTR_y)
        #np.array([[disp1, f1],[disp2, f2],[disp3, f3]]
        return utils.create_ASTR(np.array([[disp1, f1],[disp2, f2],[disp3, f3]]),negative=np.array([[-disp1, -f1],[-disp2, -f2],[-disp3, -f3]]))
    def dataframe(self):
        return pd.DataFrame([[self.s_1,self.s_2,self.s_3,self.slip2tau(self.s_1),self.slip2tau(self.s_2),self.slip2tau(self.s_3)]],columns=['s_1','s_2','s_3','t_1','t_2','t_3'])
    def curve_force(self,stop=9,num=1000,title='bond force–slip relationship',label1='MC2010',label2='ASTR', returnData = False):
        fig = plt.figure(figsize = (6,4))
        ax = fig.add_subplot(111)
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
        x=np.linspace(0, stop, num)
        y=[self.tau2force([self.slip2tau(i)])[0]/1000 for i in x]
        ax.plot(x,y,'-', linewidth=2, markersize=5,label=label1)
        astr_force=self.astr_force(stop=stop)
        ASTR_x = utils.ASTR_plot(astr_force, stop, 0, plotSide='positive', returnData = True)[:,0]
        ASTR_y = np.array(utils.ASTR_plot(astr_force, stop, 0, plotSide='positive', returnData = True)[:,1])/1000
        ax.plot(ASTR_x,ASTR_y,'-', linewidth=2, markersize=5,label=label2)
        ax.legend()
        ax.set_title(title)
        ax.set_xlabel('Slip [mm]')
        ax.set_ylabel('Force [kN]')
        ax.set_xlim(0,None)
        ax.set_ylim(0,None)
        plt.show()
        if returnData: return [x,y]
    def adaptic_print(self): # real values
        s_e=0.5*self.s_1
        return utils.str_joint([self.ID,'bond', s_e, self.s_1, self.s_2, self.s_3,self.slip2tau(s_e),self.slip2tau(self.s_1),self.slip2tau(self.s_3)])

# stmdl2 to con1 conversion
# ec0t - Ec0 # initial compressive strength
# ec0 - Ec1 # secant compressive stiffness
# muec1 - Ec2 # compressive softening secant stiffness
# et0 - Et1 # secant tensile stiffness
# muet1 - Et2 # tensile softening secant stiffness
# stresc0 - fc1 # peak compressive strength
# stresc1 - fc2 # residual compressive strength
# ft # peak tensile strength
# strnc0 - epsilon_c1 # strain at peak compressive strength
# strnc1 - epsilon_c2 # strain at residual compressive strength
# strnt0 - epsilon_t1 # strain at peak tensile strength
# strnt1 - epsilon_t2 # strain when tensile stress reaches 0 ?
# pseto # plastic strain in compression
# crkso # plastic strain in tension
# pseto # plastic strain in compression at the start of the step, represents the intersection of the unloading branch with the strain axis
# crkso # plastic strain in tension at the start of the step, represents the intersection of the unloading branch with the strain axis

class stmdl2:
    # this is con1 ADAPTIC model
    def __init__(self,ec0,muec1,strnc1,stresc1,et0,muet1,strnt1,alphac,alphat,pseto=0,crkso=0,ID='stmdl2'):#pseto,crkso,
        # This subroutine calculates the stress at a monitoring point for
        # material MODEL(2).

        # Establish the stress depending on the sign of the applied
        # strain relative to the initial plastic set and crack strain
        self.ec0=ec0 # secant compressive stiffness
        self.muec1=muec1 # compressive softening secant stiffness
        self.strnc1=strnc1 # strain at residual compressive strength
        self.stresc1=stresc1 # residual compressive strength
        self.et0=et0 # secant tensile stiffness
        self.muet1=muet1 # secant tensile softening secant stiffness
        self.strnt1=strnt1 # strain when tensile stress reaches 0
        self.pset=pseto # plastic strain in compression
        self.crks=crkso # plastic strain in tension
        self.pset=pseto # plastic strain in compression at the start of the step, represents the intersection of the unloading branch with the strain axis
        self.crks=crkso # plastic strain in tension at the start of the step, represents the intersection of the unloading branch with the strain axis
        self.alphac=alphac
        self.alphat=alphat

        self.ID = ID
        # Derived, not used in stress
        self.ec0t=(1+np.abs(alphac))*ec0 # secant compressive stiffness
        self.strnc0=(stresc1-muec1*strnc1)/(ec0-muec1) # strain at peak compressive strength
        self.stresc0=self.strnc0*self.ec0 # peak compressive strength
        self.strnt0=-muet1*strnt1/(et0-muet1) # strain at peak tensile strength
        if alphat>0: self.et0t=(1+np.abs(alphat))*et0 # secant tensile stiffness
        else: self.et0t=et0
        self.ft=self.et0*self.strnt0 # peak tensile strength

        data = np.array([[self.ID,self.stresc0, self.stresc1, self.ft, self.ec0t, self.ec0,
                          self.muec1, self.et0, self.muet1, self.strnc0,
                          self.strnc1, self.strnt0, self.strnt1, self.alphac, self.alphat]])
        self.prop = pd.DataFrame(data,index=data[:,0])
        self.prop.columns = ['$$ID$$','$$f_{c1}[MPa]$$','$$f_{c2}[MPa]$$', '$$f_{t}[MPa]$$','$$E_{c0}[MPa]$$','$$E_{c1}[MPa]$$','$$E_{c2}[MPa]$$',
            '$$E_{t1}[MPa]$$','$$E_{t2}[MPa]$$','$$e_{c1}$$','$$e_{c2}$$','$$e_{t1}$$', '$$e_{t2}$$', '$$alpha_{c}$$', '$$alpha_{t}$$']

    @classmethod # alternative constructor
    def from_ADAPTIC(cls, ec1,fc1,ec2,fc2,et1,ft,et2,alphac,alphat, ID='stmdl2'):
        strnc1=-fc1/ec1+(fc1-fc2)/ec2
        strnt1=ft/et1-ft/et2
        return cls(ec0=ec1,muec1=ec2,strnc1=strnc1,stresc1=-fc2,et0=et1,muet1=et2,strnt1=strnt1,alphac=alphac,alphat=alphat, ID = ID)

    def y_E1(self,x,E1,E0,x0,x1,y0,printing=True):
        a=(E0-E1)/2/(x0-x1)
        b=E0-2*a*x0
        c=y0-a*x0**2-b*x0
        if printing: print('y(x) = {0}x**2+{1}x+{2}'.format(a,b,c))
        return (x-x0)**2*(E0-E1)/(2*(x0-x1))+E0*(x-x0)+y0

    def y_S(self,x,S,E0,x0,x1,y0,printing=True):
        a=2*(E0-S)/2/(x0-x1)
        b=E0-2*a*x0
        c=y0-a*x0**2-b*x0
        if printing: print('y(x) = {0}x**2+{1}x+{2}'.format(a,b,c))
        return (x-x0)**2*(E0-S)/(x0-x1)+E0*(x-x0)+y0

    def y_prime_E1(self,x,E1,E0,x0,x1,y0,printing=True):
        a=(x*(E0-E1)+x0*E1-x1*E0)/(x0-x1)
        b=(x**2*(E1-E0)-x0**2*(E0+E1)+2*x0*x1*E0)/(2*(x0-x1))+y0
        if printing: print('y\'(x) = {0}x+{1}'.format(a,b))
        return (a,b)

    def y_prime_S(self,x,S,E0,x0,x1,y0,printing=True):
        a=(x*2*(E0-S)+x0*(2*S-E0)-x1*E0)/(x0-x1)
        b=(S*x**2-S*x0**2-x**2*E0+x0*x1*E0)/(x0-x1)+y0
        if printing: print('y\'(x) = {0}x+{1}'.format(a,b))
        return (a,b)

    def secant(self,S,x0,y0,printing=True):
        a=S
        b=y0-a*x0
        if printing: print('secant: y(x) = {0}x+{1}'.format(a,b))
        return (a,b)

    def stress(self,strn,printing=False,retn='stress'):
        self.strn=strn

        ec0=self.ec0 # Secant compressive stiffness
        muec1=self.muec1 # Compressive softening stiffness
        strnc1=self.strnc1 # strain at residual compressive strength
        stresc1=self.stresc1 # residual compressive strength
        et0=self.et0 # Tensile stiffness
        muet1=self.muet1 # Tensile softening stiffness
        strnt1=self.strnt1 # strain at peak tensile strength ?
        pseto=self.pset # plastic strain in compression
        crkso=self.crks # plastic strain in tension
        alphac=self.alphac
        alphat=self.alphat

        # Establish the stress depending on the sign of the applied
        # strain relative to the initial plastic set and crack strain

        if(strn==pseto):
            stres=0.0
            pset=pseto
            crks=crkso
            etan=0

        if(strn<=pseto):     # Compressive strain increment

            # NOTE: alphac is the relative difference between the initial
            #       compressive tangent modulus and the secant modulus

            ec0t=(1+np.abs(alphac))*ec0
            if(alphac!=0): strnc0=(stresc1-muec1*strnc1)/(ec0-muec1)

            # Obtain the stress assuming elastic conditions, and determine
            # the force from the limiting curve

            strese=ec0t*(strn-pseto)

            if(alphac!=0 and strn>strnc0):
                stresl=strn*(ec0t+(ec0-ec0t)*strn/strnc0)
                etan=ec0t+2*(ec0-ec0t)*strn/strnc0

            elif(strn>strnc1):
                dstrn1=strn-strnc1

                # linear softening branch
                stresl=stresc1+muec1*dstrn1
                etan=muec1

                if(alphac<0):
                    # overaly with cubic function
                    dstrn0=strnc0-strnc1
                    dstrn1=dstrn1/dstrn0
                    stresl=stresl+2*alphac*muec1*dstrn0*dstrn1*(dstrn1-1)*(dstrn1-0.5)
                    etan=etan+alphac*muec1*(1-6*dstrn1*(1-dstrn1))

            else:
                # residual compressive strength
                stresl=stresc1
                etan=0.0

            # Establish the stress and the plastic set

            if(strese>stresl):
                stres=strese
                pset=pseto
                etan=ec0t

            else:
                stres=stresl
                pset=strn-stresl/ec0t

            crks=crkso

        elif(et0==0.0 or strn>=pseto+strnt1 or crkso>=strnt1):   # No tensile resistance
            stres=0.0
            pset=pseto
            crks=strnt1
            etan=0

        else:                                         # Tensile strain increment
            # NOTE: The tensile response is modified so that unloading points
            #       towards the origin of compressive loading response (i.e.
            #       plastic compressive strain), and the cracking strain is
            #       now defined as the maximum strain relative to the origin
            #       (rather than the unloading strain)

            pset=pseto

            strnt0=-muet1*strnt1/(et0-muet1)

            # Obatin relevant tensile strain to establish current stress on
            # loading/softening envelope

            dstrn=strn-pseto
            onEnv=True

            if(dstrn<=crkso): # unloading curve
                dstrn=crkso
                onEnv=False      # Use secant stiffness

            # Obtain relevant stress and tangent modulus on envelope

            if(dstrn<=strnt0):   # loading envelope
                if(alphat<=0):
                    # linear envelope: simple linear case, no need for further checks
                    stres=et0*dstrn # before et0*(strn-pseto)
                    stres=et0*(strn-pseto)
                    etan=et0
                    crks=crkso

                else:
                    # quadratic envelope
                    et0t=(1+alphat)*et0
                    stres=dstrn*(et0t+(et0-et0t)*dstrn/strnt0)
                    if(onEnv): etan=et0t+2*(et0-et0t)*dstrn/strnt0

            else:                     # softening envelope

                dstrn1=dstrn-strnt1

                # linear softening branch
                stres=muet1*dstrn1
                if(onEnv): etan=muet1

                if(alphat!=0):

                    dstrn0=strnt0-strnt1
                    dstrn1=dstrn1/dstrn0

                    if(alphat>0):
                        # overlay linear with cubic function
                        stres=stres-2*alphat*muet1*dstrn0*dstrn1*(dstrn1-1)*(dstrn1-0.5)
                        if(onEnv): etan=etan-alphat*muet1*(1-6*dstrn1*(1-dstrn1))
                    else:
                        # overlay linear with quadratic function
                        stres=stres-alphat*muet1*dstrn0*dstrn1*(dstrn1-1)
                        if(onEnv): etan=etan-alphat*muet1*(2*dstrn1-1)

            crks=dstrn

            # Determine tangent modulus as secant unloading stiffness to
            # compression origin if stress state is not on envelope

            if(not onEnv):

                etan=stres/dstrn # because dstrn=crkso
                stres=etan*(strn-pseto)

        self.strn=strn
        self.stres=stres
        self.etan=etan
        self.crks=crks
        self.pset=pset

        if retn=='etan': return self.etan
        elif retn=='all': return self.stres,self.etan
        else: return self.stres

    def plot(self,strain='',retn='stress',title='con1',lineType='-',legend=True,lbl='stmdl2',xlim=(None,None),ylim=(None,None),ylabel='Stress [MPa]',xlabel='Strain',pseto='',crkso=''):
        if strain=='':strain=np.arange(1.05*np.absolute(self.strnc1),1.05*np.absolute(self.strnt1),0.00001)
        fig = plt.figure(figsize = (6,4))
        ax = fig.add_subplot(111)
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
        stress,etan=[],[]
        for j,i in enumerate(strain):
            if crkso!='': self.crks=crkso
            if pseto!='': self.pset=pseto
            self.stress(i)
            stress.append(self.stres)
            etan.append(self.etan)
#             if j>2:
#                 X=strain[-3:-1]
#                 Y=stress[-3:-1]
#                 print('step: {0}, etan: {1}, slope_intercept: {2}'.format(j,self.etan,self.slope_intercept(X,Y)[0]))
        stress=np.array(stress)
        etan=np.array(etan)
        strain=strain.reshape(len(strain),1)
        stress=stress.reshape(len(stress),1)
        etan=etan.reshape(len(etan),1)
        self.df=pd.DataFrame(np.hstack((stress,strain,etan)),columns=['stress','strain','etan'])
        if 'etan' in retn:
            if retn=='etan' or retn=='etan1':ax.plot(self.df['strain'],self.df['etan'],lineType, linewidth=2, markersize=5,label=lbl)
            if retn=='etan' or retn=='etan2':
                slope=self.slope_intercept(strain,stress)
                ax.plot(strain[:-1],slope,lineType, linewidth=2, markersize=5,label='slope')
        else: ax.plot(self.df['strain'],self.df['stress'],lineType, linewidth=2, markersize=5,label=lbl)
        if legend: ax.legend()
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.show()

    def slope_intercept(self,X,Y):
        slope=[]
        for i in range(0,len(X)-1):
            if (X[i+1] - X[i])==0:
                print('ERROR: slope_intercept({},{})'.format([X[i+1],X[i]],[Y[i+1],Y[i]]))
                slope.append(math.nan)
            else: slope.append((Y[i+1] - Y[i]) / (X[i+1] - X[i]) )
        return slope
