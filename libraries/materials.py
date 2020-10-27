import numpy as np
import pandas as pd
import utils
import matplotlib.pyplot as plt

class con1:
    def __init__(self, ID, fc1, length, epsilon_t2 = 0.001, fc2_factor = 0.1, ft_factor = 1, characteristic = True, Ec2 = '',Et2 = '',strain_prec=5):
        self.resid_str = fc2_factor
        self.ID = ID
        self.fc1 = round(fc1, 1)
        self.fc2 = round(self.resid_str * fc1, 1)
        self.length = length
        if characteristic:
            self.fcm = round(fc1+8,2)
        else:
            self.fcm = fc1
        if fc1 <= 50:
            self.ft = round(ft_factor * 0.3 * self.fcm ** (2/3), 1)
        else:
            self.ft = round(ft_factor * 2.12*np.log(1+0.1*self.fcm), 1)
        self.Gf = round(73 * self.fcm**0.18/1000, 3)
        self.Ec0 = round(int(21500*(self.fcm/10)**(1/3)),-2)
        self.poisson = 0.2
        self.Gc = round(250 * self.Gf, 1)
        self.epsilon_1c = round(5 * fc1 / self.Ec0 /3, strain_prec)
        self.Ec1 = int(round(fc1 / self.epsilon_1c, -2))
        self.alpha = min(max(0,round((self.Ec0 - self.Ec1)/self.Ec1,2)),1)
        if Ec2 != '':
            self.Ec2 = Ec2
            self.epsilon_2c = round((self.fc1-self.fc2)/-Ec2+self.epsilon_1c, strain_prec)
        else:
            self.epsilon_2c = round(self.epsilon_1c + 3 * self.Gc / (2 * length * fc1), strain_prec)
            self.Ec2 = - int(round((1 - self.resid_str) * fc1 /(self.epsilon_2c - self.epsilon_1c), -2))
        self.Et1 = self.Ec0
        self.epsilon_1t = round(self.ft / self.Et1, strain_prec)
        if Et2 != '':
            self.Et2 = Et2
            self.epsilon_2t = round((self.ft)/-Et2+self.epsilon_1t, strain_prec)
        else:
            self.epsilon_2t = epsilon_t2
            self.Et2 = - int(round(self.ft /(self.epsilon_2t - self.epsilon_1t), -2))

    def adaptic_print(self):
        line = utils.str_joint([self.ID,'con1', self.Ec1, self.fc1, self.Ec2, self.fc2, self.Et1,
                                      self.ft, self.Et2, self.alpha])
        return line

    def data_frame(self):
        data = np.array([[self.ID, self.length, self.fc1, self.fc2, self.ft, self.Ec0, self.Ec1,
                          self.Ec2, self.Et1, self.Et2, self.Gf, self.Gc, self.epsilon_1c,
                          self.epsilon_2c, self.epsilon_1t, self.epsilon_2t,  self.alpha]])
        df = pd.DataFrame(data,index=data[:,0])
        df.columns = ['ID', '$$h[mm]$$', '$$f_{c1}[MPa]$$','$$f_{c2}[MPa]$$', '$$f_{t}[MPa]$$',
                      '$$E_{c0}[MPa]$$','$$E_{c1}[MPa]$$','$$E_{c2}[MPa]$$','$$E_{t1}[MPa]$$',
                      '$$E_{t2}[MPa]$$','$$G_{f}[N/mm]$$','$$G_{c}[N/mm]$$','$$e_{c1}$$',
                      '$$e_{c2}$$','$$e_{t1}$$', '$$e_{t2}$$', '$$alpha$$']
        return df


class stl1:
    def __init__(self, ID, E1, fy, fu, epsilon_u):
        self.ID = ID
        self.E1 = E1
        self.fy = fy
        self.fu = fu
        self.epsilon_u = round(epsilon_u,3)
        self.epsilon_y = round(fy / E1,4)
        self.E2 = round((fu - fy) / (epsilon_u - self.epsilon_y),1)
        self.mu = round(self.E2 / E1,7)

    def adaptic_print(self):
        line = utils.str_joint([self.ID,'stl1', self.E1, self.fy, self.mu])
        return line

    def data_frame(self):
        data = np.array([[self.ID, self.E1, self.E2, self.fy, self.fu, self.epsilon_y, self.epsilon_u,
                          self.mu]])
        df = pd.DataFrame(data,index=data[:,0])
        df.columns = ['ID', '$$E_{1}[MPa]$$', '$$E_{2}[MPa]$$', '$$f_{y}[MPa]$$', '$$f_{u}[MPa]$$',
                      '$$e_{y}$$','$$e_{u}$$','$$mu$$']
        return df

    def stress_df(self):
        df = pd.DataFrame([[0,0],[self.epsilon_y,self.fy],[self.epsilon_u,self.fu]],columns=['strain','stress'])
        return df

class bond:
    def __init__(self,c,f_cm,L,dia,n_bars,case=1,redFact=1):
        # case - refer to Table 6.1-1 MC2010:
        # 1 - Pull-out, good bond
        # 2 - Pull-out, all other bond cond
        # 3 - Splitting, good bond cond, unconfined
        # 4 - Splitting, good bond cond, stirrups
        # 5 - Splitting, all other bond cond, unconfined
        # 6 - Splitting, all other bond cond, stirrups
        self.case=case
        self.f_cm = f_cm # MPa
        self.L=L
        self.dia=dia
        self.n_bars=n_bars
        self.redFact=redFact
        if self.case ==1:
            self.tau_max = 2.5 * self.f_cm**0.5 * redFact
            self.s_1 = 1 # mm
            self.s_2 = 2 # mm
            self.s_3 = c # mm, clear distance between ribs
            self.alpha = 0.4
            self.tau_bf = 0.4 * self.tau_max
        elif self.case ==4:
            self.tau_max = 2.5 * self.f_cm**0.5
            self.tau_bu_split=8*(f_cm/25)**0.25 * redFact
            self.s_1 = 1/self.tau_max*self.tau_bu_split # mm
            #self.s_2 = self.s_1 # mm
            self.s_3 = 0.5 * c # mm
            self.alpha = 0.4
            self.tau_bf = 0.4 * self.tau_bu_split
        else: print('Case error')
    def slip2tau(self,s):
        if self.case ==1:
            if 0 <= s <= self.s_1:
                tau = self.tau_max*(s/self.s_1)**self.alpha
            elif self.s_1 <= s <= self.s_2:
                tau = self.tau_max
            elif self.s_2 <= s <= self.s_3:
                tau = self.tau_max - (self.tau_max-self.tau_bf)*(s-self.s_2)/(self.s_3-self.s_2)
            else:
                tau = self.tau_bf
        elif self.case ==4:
            if 0 <= s <= self.s_1:
                tau = self.tau_bu_split*(s/self.s_1)**self.alpha
            elif self.s_1 <= s <= self.s_3:
                tau = self.tau_bu_split - (self.tau_bu_split-self.tau_bf)*(s-self.s_1)/(self.s_3-self.s_1)
            else:
                tau = self.tau_bf
        return tau
    def force2tau(self,force):
        return force/(self.n_bars*self.dia*np.pi*self.L)
    def tau2force(self,tau):
        return tau*(self.n_bars*self.dia*np.pi*self.L)
    def curve(self,stop=9,num=50,title='bond stress–slip relationship'):
        fig = plt.figure(figsize = (6,4))
        ax = fig.add_subplot(111)
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
        x=np.linspace(0, stop, num)
        y=[self.slip2tau(i) for i in x]
        ax.plot(x,y,'-', linewidth=2, markersize=5)
        ax.set_title(title)
        ax.set_xlabel('Slip [mm]')
        ax.set_ylabel('Stress [MPa]')
        ax.set_xlim(0,None)
        ax.set_ylim(0,None)
        plt.show()
    def astr_curve(self):
        disp1 = self.s_1
        disp2 = self.s_3
        disp3 = self.s_3+1
        f1 = self.tau2force(self.slip2tau(disp1))
        f2 = self.tau2force(self.slip2tau(disp2))
        f3 = self.tau2force(self.slip2tau(disp3))
        #np.array([[disp1, f1],[disp2, f2],[disp3, f3]]
        return f1,f2,f3,disp1,disp2,disp3
    def dataframe(self):
        return pd.DataFrame([[self.s_1,self.s_2,self.s_3,self.slip2tau(self.s_1),self.slip2tau(self.s_2),self.slip2tau(self.s_3)]],columns=['s_1','s_2','s_3','t_1','t_2','t_3'])
    def curve_force(self,stop=9,num=50,title='bond force–slip relationship'):
        fig = plt.figure(figsize = (6,4))
        ax = fig.add_subplot(111)
        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')
        x=np.linspace(0, stop, num)
        y=[(self.tau2force(self.slip2tau(i)))/1000 for i in x]
        ax.plot(x,y,'-', linewidth=2, markersize=5,label='MC2010')
        f1,f2,f3,disp1,disp2,disp3=self.astr_curve()
        ax.plot([0,disp1,disp2,disp3],[0,f1/1000,f2/1000,f3/1000],'-', linewidth=2, markersize=5,label='ASTR')
        ax.legend()
        ax.set_title(title)
        ax.set_xlabel('Slip [mm]')
        ax.set_ylabel('Force [kN]')
        ax.set_xlim(0,None)
        ax.set_ylim(0,None)
        plt.show()
