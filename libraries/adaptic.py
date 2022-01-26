import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# ADAPTIC classes
# =============================================================================
class adap0: # base class with funtions to read num files
    def __init__(self,name):
        self.name = name

    def readFile(cls, title, cutoff = None,folderPath="",numPath=""): # numPath = str('/num/'), path='data/'
        # wordsNum=[]
        # for data in open(folderPath+numPath+title+str(".num"),'r'):
        #     wordsNum.append(data.split())
        # if cutoff:
        #     for i in range(len(wordsNum)):
        #         try:
        #             if wordsNum[i][0] == '#io1' and wordsNum[i+2][0] == str(cutoff):
        #                 break
        #         except:
        #             pass
        #     wordsNum = wordsNum[0:i-2]
        # return wordsNum

        #update 2021/07/26, cutoff specific steps
        # wordsNum=[]
        # for data in open(folderPath+numPath+title+str(".num"),'r'):
        #     wordsNum.append(data.split())
        # if cutoff:
        #     lineNo = []
        #     for i in range(len(wordsNum)):
        #         try:
        #             if wordsNum[i][0] == '#io1':
        #                 lineNo.append([int(wordsNum[i+2][0]),i])
        #         except:
        #             pass
        #     cutoff.sort(reverse=True)
        #     for i in cutoff:
        #         if i==lineNo[-1][0]:
        #             wordsNum=wordsNum[:(lineNo[-1][1]-2)]
        #         elif i<lineNo[-1][0]:
        #             del wordsNum[(lineNo[i-1][1]):(lineNo[i][1])]
        # return wordsNum

        # update 2021/11/17, cutoff specific steps
        wordsNum=[]
        for data in open(folderPath+numPath+title+str(".num"),'r'):
            wordsNum.append(data.split())
        if cutoff:
            lineNo = []
            for i in range(len(wordsNum)):
                try:
                    if wordsNum[i][0] == '#io1':
                        lineNo.append([int(wordsNum[i+2][0]),i])
                except:
                    pass
            steps=list(list(zip(*lineNo))[0])
            line=list(list(zip(*lineNo))[1])
            if cutoff<=0: del wordsNum[(line[steps[cutoff]-1]-2):]
            elif steps[-1] > cutoff:
                del wordsNum[(line[steps.index(cutoff)]-2):]
        return wordsNum

    def convertNum(self, wordsNum, phrase, startRow, column, convType):
        variable=[]
        iRow=-1
        for i in wordsNum:
            iRow += 1
            if i==[phrase]:
                tempVector=[]
                jRow = iRow + startRow
                while wordsNum[jRow] != []:
                    if convType == "int":
                        tempVal = int(wordsNum[jRow][column])
                    elif convType == "float":
                        try: tempVal = float(wordsNum[jRow][column])
                        except: tempVal = float(0)
                    elif convType == "str":
                        tempVal = wordsNum[jRow][column]
                    tempVector.append(tempVal)
                    jRow += 1
                tempVector=np.array(tempVector)
                if len(variable)==0:
                    variable = tempVector
                else:
                    variable=np.column_stack((variable,tempVector))
                if convType == "str":
                    break
        if convType != "str":
            variable=np.column_stack((np.zeros((variable.shape[0],1)),variable))
        return variable

class adap1(adap0): # base adaptic class processing the num
    def __init__(self,title, cutoff, folderPath, numPath):
        super(adap1, self).__init__(title)
        self.wordsNum = self.readFile(title, cutoff,folderPath,numPath)
        self.step = self.convertNum(self.wordsNum,"#io1",2,0,"int")
        self.nodeName = self.convertNum(self.wordsNum,"#in2",3,0,"str")

    def returnVec(self,attribute,idx=':'):
        self.hasAttribute(attribute)
        return eval('self.'+attribute+'['+str(idx)+'].T')

    # def pseudoStatic(self, LF, dispMax):
    #     pseudoCrv = [0]
    #     for i in range(1,LF.shape[1]):
    #         pseudoTemp=np.trapz(LF[0,0:i+1],dispMax[0,0:i+1])/dispMax[0,i]
    #         pseudoCrv.append(pseudoTemp)
    #     return np.array(pseudoCrv)[None, :]

    def restrainedName_create(self):
        self.restrainedName = self.convertNum(self.wordsNum,"#in1",3,0,"str")

    def cbpName_create(self):
        self.cbpName = self.convertNum(self.wordsNum,"#ie1",3,0,"str")

    def jelName_create(self):
        self.jelName = self.convertNum(self.wordsNum,"#ie11",3,0,"str")

    def lnkName_create(self):
        self.lnkName = self.convertNum(self.wordsNum,"#ie18",3,0,"str")

    def nodeDispX_create(self):
        self.nodeDispX = self.convertNum(self.wordsNum,"#in2",3,1,"float")

    def nodeDispY_create(self):
        self.nodeDispY = self.convertNum(self.wordsNum,"#in2",3,2,"float")

    def nodeDispY_create(self):
        self.nodeDispY = self.convertNum(self.wordsNum,"#in2",3,2,"float")

    def restrainedX_create(self):
        self.hasAttribute("restrainedName")
        self.restrainedX = self.convertNum(self.wordsNum,"#in1",3,1,"float")

    def restrainedY_create(self):
        self.hasAttribute("restrainedName")
        self.restrainedY = self.convertNum(self.wordsNum,"#in1",3,2,"float")

    # def time_create(self):
    #     self.time = self.convertNum(self.wordsNum,"#io1",2,2,"float")

    # def LF_create(self):
    #     self.LF = self.convertNum(self.wordsNum,"#io1",2,2,"float")

    def hasAttribute(self, att):
        if hasattr(self, att):
            pass
        else:
            getattr( self, str(att+"_create"))()

class adaptic2D(adap1):
    def __init__(self,title, cutoff = None,folderPath="",numPath=""):
        super(adaptic2D, self).__init__(title, cutoff, folderPath, numPath)

    def gaussName_create(self):
        temp1 = self.convertNum(self.wordsNum,"#ie1s",3,0,"str")
        temp2 = self.convertNum(self.wordsNum,"#ie1s",3,1,"str")
        self.gaussName = np.array([temp1[i]+"_"+temp2[i] for i in range(len(temp2))])
        del temp1, temp2

    def nodeDispRZ_create(self):
        self.nodeDispRZ = self.convertNum(self.wordsNum,"#in2",3,3,"float")

    def cbpM1_create(self):
        self.hasAttribute("cbpName")
        self.cbpM1 = self.convertNum(self.wordsNum,"#ie1",3,1,"float")

    def cbpTheta1_create(self):
        self.hasAttribute("cbpName")
        self.cbpTheta1 = self.convertNum(self.wordsNum,"#ie1d1",3,1,"float")

    def restrainedRZ_create(self):
        self.hasAttribute("restrainedName")
        self.restrainedRZ = self.convertNum(self.wordsNum,"#in1",3,3,"float")

    def cbpM2_create(self):
        self.hasAttribute("cbpName")
        self.cbpM2 = self.convertNum(self.wordsNum,"#ie1",3,2,"float")

    def cbpTheta2_create(self):
        self.hasAttribute("cbpName")
        self.cbpTheta2 = self.convertNum(self.wordsNum,"#ie1d1",3,2,"float")

    def cbpF_create(self):
        self.hasAttribute("cbpName")
        self.cbpF = self.convertNum(self.wordsNum,"#ie1",3,3,"float")

    def cbpDelta_create(self):
        self.hasAttribute("cbpName")
        self.cbpDelta = self.convertNum(self.wordsNum,"#ie1d1",3,3,"float")

    def jelF_create(self):
        self.hasAttribute("jelName")
        self.jelF = self.convertNum(self.wordsNum,"#ie11",3,1,"float")

    def jelV_create(self):
        self.hasAttribute("jelName")
        self.jelV = self.convertNum(self.wordsNum,"#ie11",3,2,"float")

    def jelM_create(self):
        self.hasAttribute("jelName")
        self.jelM = self.convertNum(self.wordsNum,"#ie11",3,3,"float")

    def lnkM1_create(self):
        self.hasAttribute("lnkName")
        self.lnkM1 = self.convertNum(self.wordsNum,"#ie18",3,1,"float")

    def lnkM2_create(self):
        self.hasAttribute("lnkName")
        self.lnkM2 = self.convertNum(self.wordsNum,"#ie18",3,2,"float")

    def lnkF_create(self):
        self.hasAttribute("lnkName")
        self.lnkF = self.convertNum(self.wordsNum,"#ie18",3,3,"float")

    def gauss1StrainB_create(self):
        self.hasAttribute("gaussName")
        self.gauss1StrainB = self.convertNum(self.wordsNum,"#ie1s",3,2,"float")

    def gauss1StrainT_create(self):
        self.hasAttribute("gaussName")
        self.gauss1StrainT = self.convertNum(self.wordsNum,"#ie1s",3,4,"float")

    def gauss1StrainAv_create(self):
        self.hasAttribute("gauss1StrainT")
        self.hasAttribute("gauss1StrainB")
        self.hasAttribute("gaussName")
        self.gauss1StrainAv = (self.gauss1StrainT+self.gauss1StrainB)/2

    def gauss2StrainAv_create(self):
        self.hasAttribute("gauss2StrainT")
        self.hasAttribute("gauss2StrainB")
        self.hasAttribute("gaussName")
        self.gauss2StrainAv = (self.gauss2StrainT+self.gauss2StrainB)/2

    def gauss1StressAv_create(self):
        self.hasAttribute("gauss1StressT")
        self.hasAttribute("gauss1StressB")
        self.hasAttribute("gaussName")
        self.gauss1StressAv = (self.gauss1StressT+self.gauss1StressB)/2

    def gauss2StressAv_create(self):
        self.hasAttribute("gauss2StressT")
        self.hasAttribute("gauss2StressB")
        self.hasAttribute("gaussName")
        self.gauss2StressAv = (self.gauss2StressT+self.gauss2StressB)/2

    def gauss1StressB_create(self):
        self.hasAttribute("gaussName")
        self.gauss1StressB = self.convertNum(self.wordsNum,"#ie1s",3,3,"float")

    def gauss1StressT_create(self):
        self.hasAttribute("gaussName")
        self.gauss1StressT = self.convertNum(self.wordsNum,"#ie1s",3,5,"float")

    def gauss2StrainB_create(self):
        self.hasAttribute("gaussName")
        self.gauss2StrainB = self.convertNum(self.wordsNum,"#ie1s",3,6,"float")

    def gauss2StrainT_create(self):
        self.hasAttribute("gaussName")
        self.gauss2StrainT = self.convertNum(self.wordsNum,"#ie1s",3,8,"float")

    def gauss2StressB_create(self):
        self.hasAttribute("gaussName")
        self.gauss2StressB = self.convertNum(self.wordsNum,"#ie1s",3,7,"float")

    def gauss2StressT_create(self):
        self.hasAttribute("gaussName")
        self.gauss2StressT = self.convertNum(self.wordsNum,"#ie1s",3,9,"float")

    def bndName_create(self):
        self.bndName = self.convertNum(self.wordsNum,"#ie23",3,0,"str")

    def bndM1_create(self):
        self.hasAttribute("bndName")
        self.bndM1 = self.convertNum(self.wordsNum,"#ie23",3,1,"float")

    def bndM2_create(self):
        self.hasAttribute("bndName")
        self.bndM2 = self.convertNum(self.wordsNum,"#ie23",3,2,"float")

    def bndF_create(self):
        self.hasAttribute("bndName")
        self.bndF = self.convertNum(self.wordsNum,"#ie23",3,3,"float")

    def bndGaussName_create(self):
        temp1 = self.convertNum(self.wordsNum,"#ie23s",3,0,"str")
        temp2 = self.convertNum(self.wordsNum,"#ie23s",3,1,"str")
        self.bndGaussName = np.array([temp1[i]+"_"+temp2[i] for i in range(len(temp2))])
        del temp1, temp2

    def bndStrain_create(self):
        self.hasAttribute("bndGaussName")
        self.bndStrain = self.convertNum(self.wordsNum,"#ie23s",3,2,"float")

    def bndStress_create(self):
        self.hasAttribute("bndGaussName")
        self.bndStress = self.convertNum(self.wordsNum,"#ie23s",3,3,"float")

    def bndSlip_create(self):
        self.hasAttribute("bndGaussName")
        self.bndSlip = self.convertNum(self.wordsNum,"#ie23s",3,4,"float")

    def bndBond_create(self):
        self.hasAttribute("bndGaussName")
        self.bndBond = self.convertNum(self.wordsNum,"#ie23s",3,5,"float")

    def cncName_create(self):
        self.cncName = self.convertNum(self.wordsNum,"#ie24",3,0,"str")

    def cncM1_create(self):
        self.hasAttribute("cncName")
        self.cncM1 = self.convertNum(self.wordsNum,"#ie24",3,1,"float")

    def cncM2_create(self):
        self.hasAttribute("cncName")
        self.cncM2 = self.convertNum(self.wordsNum,"#ie24",3,2,"float")

    def cncF_create(self):
        self.hasAttribute("cncName")
        self.cncF = self.convertNum(self.wordsNum,"#ie24",3,3,"float")

    def cncGaussName_create(self):
        temp1 = self.convertNum(self.wordsNum,"#ie24s",3,0,"str")
        temp2 = self.convertNum(self.wordsNum,"#ie24s",3,1,"str")
        self.cncGaussName = np.array([temp1[i]+"_"+temp2[i] for i in range(len(temp2))])
        del temp1, temp2

    def cncStrain_create(self):
        self.hasAttribute("cncGaussName")
        self.cncStrain = self.convertNum(self.wordsNum,"#ie24s",3,2,"float")

    def cncStress_create(self):
        self.hasAttribute("cncGaussName")
        self.cncStress = self.convertNum(self.wordsNum,"#ie24s",3,3,"float")

    def cncGamma_create(self):
        self.hasAttribute("cncGaussName")
        self.cncGamma = self.convertNum(self.wordsNum,"#ie24s",3,4,"float")

    def cncTau_create(self):
        self.hasAttribute("cncGaussName")
        self.cncTau = self.convertNum(self.wordsNum,"#ie24s",3,5,"float")

    def findIndice(self, att, ID):
        indice = "error"
        if "restrained" in att:
            att = "restrainedName"
        if "cbp" in att:
            att = "cbpName"
        if "jel" in att:
            att = "jelName"
        if "lnk" in att:
            att = "lnkName"
        if "gauss" in att:
            att = "gaussName"
        if att == "bndM1" or att == "bndM2" or att == "bndF":
            att = "bndName"
        if att == "bndStrain" or att == "bndStress" or att == "bndSlip" or att == "bndBond":
            att = "bndGaussName"
        if att == "cncM1" or att == "cncM2" or att == "cncF":
            att = "cncName"
        if att == "cncStrain" or att == "cncStress" or att == "cncGamma" or att == "cncTau":
            att = "cncGaussName"
        #print('att = {0}'.format(att))
        for i, j in enumerate(getattr(self, att)):
            if j == ID:
                indice = i
                break
        if indice == "error":
            print('indice {0} in attribute {1} does not exist'.format(ID, att))
        return indice

class adaptic3D(adap1):
    def __init__(self,title, cutoff = None,folderPath="",numPath=""):
        super(adaptic3D, self).__init__(title, cutoff, folderPath, numPath)

    def gaussName_create(self): # element - material - position
        temp1 = self.convertNum(self.wordsNum,"#ie31s",3,0,"str")
        temp2 = self.convertNum(self.wordsNum,"#ie31s",3,1,"str")
        temp2 = np.array([i.replace(".", "_") for i in temp2])
        self.gaussName = np.array([temp1[i]+"_"+temp2[i] for i in range(len(temp2))])
        del temp1, temp2

    def cvsName_create(self): # element - material - position
        temp1 = self.convertNum(self.wordsNum,"#ie52s1",3,0,"str")
        temp2 = self.convertNum(self.wordsNum,"#ie52s1",3,1,"str")
        temp2 = np.array([i.replace(".", "_") for i in temp2])
        self.cvsName = np.array([temp1[i]+"_"+temp2[i] for i in range(len(temp2))])
        del temp1, temp2

    def nodeDispZ_create(self):
        self.nodeDispZ = self.convertNum(self.wordsNum,"#in2",3,3,"float")

    def nodeDispRX_create(self):
        self.nodeDispRX = self.convertNum(self.wordsNum,"#in2",3,4,"float")

    def nodeDispRY_create(self):
        self.nodeDispRY = self.convertNum(self.wordsNum,"#in2",3,5,"float")

    def nodeDispRZ_create(self):
        self.nodeDispRZ = self.convertNum(self.wordsNum,"#in2",3,6,"float")

    def cbpMy1_create(self):
        self.hasAttribute("cbpName")
        self.cbpMy1 = self.convertNum(self.wordsNum,"#ie31",3,1,"float")

    def cbpMz1_create(self):
        self.hasAttribute("cbpName")
        self.cbpMz1 = self.convertNum(self.wordsNum,"#ie31",3,2,"float")

    def cbpMy2_create(self):
        self.hasAttribute("cbpName")
        self.cbpMy2 = self.convertNum(self.wordsNum,"#ie31",3,3,"float")

    def cbpMz2_create(self):
        self.hasAttribute("cbpName")
        self.cbpMz2 = self.convertNum(self.wordsNum,"#ie31",3,4,"float")

    def cbpF_create(self):
        self.hasAttribute("cbpName")
        self.cbpF = self.convertNum(self.wordsNum,"#ie31",3,5,"float")

    def cbpMT_create(self):
        self.hasAttribute("cbpName")
        self.cbpMT = self.convertNum(self.wordsNum,"#ie31",3,6,"float")

    def cbpThy1_create(self):
        self.hasAttribute("cbpName")
        self.cbpThy1 = self.convertNum(self.wordsNum,"#ie31d1",3,1,"float")

    def cbpThz1_create(self):
        self.hasAttribute("cbpName")
        self.cbpThz1 = self.convertNum(self.wordsNum,"#ie31d1",3,2,"float")

    def cbpThy2_create(self):
        self.hasAttribute("cbpName")
        self.cbpThy2 = self.convertNum(self.wordsNum,"#ie31d1",3,3,"float")

    def cbpThz2_create(self):
        self.hasAttribute("cbpName")
        self.cbpThz2 = self.convertNum(self.wordsNum,"#ie31d1",3,4,"float")

    def cbpDelta_create(self):
        self.hasAttribute("cbpName")
        self.cbpDelta = self.convertNum(self.wordsNum,"#ie31d1",3,5,"float")

    def cbpTheta_create(self):
        self.hasAttribute("cbpName")
        self.cbpTheta = self.convertNum(self.wordsNum,"#ie31d1",3,6,"float")

    def restrainedRZ_create(self):
        self.hasAttribute("restrainedName")
        self.restrainedRZ = self.convertNum(self.wordsNum,"#in1",3,3,"float")

    def jelFx_create(self):
        self.hasAttribute("jelName")
        self.jelFx = self.convertNum(self.wordsNum,"#ie41",3,1,"float")

    def jelFy_create(self):
        self.hasAttribute("jelName")
        self.jelFy = self.convertNum(self.wordsNum,"#ie41",3,2,"float")

    def jelFz_create(self):
        self.hasAttribute("jelName")
        self.jelFz = self.convertNum(self.wordsNum,"#ie41",3,3,"float")

    def jelMx_create(self):
        self.hasAttribute("jelName")
        self.jelMx = self.convertNum(self.wordsNum,"#ie41",3,4,"float")

    def jelMy_create(self):
        self.hasAttribute("jelName")
        self.jelMy = self.convertNum(self.wordsNum,"#ie41",3,5,"float")

    def jelMz_create(self):
        self.hasAttribute("jelName")
        self.jelMz = self.convertNum(self.wordsNum,"#ie41",3,6,"float")

    def lnkMy1_create(self):
        self.hasAttribute("lnkName")
        self.lnkMy1 = self.convertNum(self.wordsNum,"#ie18",3,1,"float")

    def lnkMz1_create(self):
        self.hasAttribute("lnkName")
        self.lnkMz1 = self.convertNum(self.wordsNum,"#ie18",3,2,"float")

    def lnkMy2_create(self):
        self.hasAttribute("lnkName")
        self.lnkMy2 = self.convertNum(self.wordsNum,"#ie18",3,3,"float")

    def lnkMz2_create(self):
        self.hasAttribute("lnkName")
        self.lnkMz2 = self.convertNum(self.wordsNum,"#ie18",3,4,"float")

    def lnkF_create(self):
        self.hasAttribute("lnkName")
        self.lnkF = self.convertNum(self.wordsNum,"#ie18",3,5,"float")

    def lnkMT_create(self):
        self.hasAttribute("lnkName")
        self.lnkMT = self.convertNum(self.wordsNum,"#ie18",3,6,"float")

    def gauss1Strain_create(self):
        self.hasAttribute("gaussName")
        self.gauss1Strain = self.convertNum(self.wordsNum,"#ie31s",3,2,"float")

    def gauss2Strain_create(self):
        self.hasAttribute("gaussName")
        self.gauss2Strain = self.convertNum(self.wordsNum,"#ie31s",3,4,"float")

    def gauss1Stress_create(self):
        self.hasAttribute("gaussName")
        self.gauss1Stress = self.convertNum(self.wordsNum,"#ie31s",3,3,"float")

    def gauss2Stress_create(self):
        self.hasAttribute("gaussName")
        self.gauss2Stress = self.convertNum(self.wordsNum,"#ie31s",3,5,"float")

    def cvsNx_create(self):
        self.hasAttribute("cvsName")
        self.cvsNx = self.convertNum(self.wordsNum,"#ie52s1",3,2,"float")

    def cvsNy_create(self):
        self.hasAttribute("cvsName")
        self.cvsNy = self.convertNum(self.wordsNum,"#ie52s1",3,3,"float")

    def cvsNxy_create(self):
        self.hasAttribute("cvsName")
        self.cvsNxy = self.convertNum(self.wordsNum,"#ie52s1",3,4,"float")

    def cvsMx_create(self):
        self.hasAttribute("cvsName")
        self.cvsMx = self.convertNum(self.wordsNum,"#ie52s1",3,5,"float")

    def cvsMy_create(self):
        self.hasAttribute("cvsName")
        self.cvsMy = self.convertNum(self.wordsNum,"#ie52s1",3,6,"float")

    def cvsMxy_create(self):
        self.hasAttribute("cvsName")
        self.cvsMxy = self.convertNum(self.wordsNum,"#ie52s1",3,7,"float")

    def cvsQxz_create(self):
        self.hasAttribute("cvsName")
        self.cvsQxz = self.convertNum(self.wordsNum,"#ie52s1",3,8,"float")

    def cvsQyz_create(self):
        self.hasAttribute("cvsName")
        self.cvsQyz = self.convertNum(self.wordsNum,"#ie52s1",3,9,"float")

    def cvsEpsx_create(self):
        self.hasAttribute("cvsName")
        self.cvsEpsx = self.convertNum(self.wordsNum,"#ie52s2",3,2,"float")

    def cvsEpsy_create(self):
        self.hasAttribute("cvsName")
        self.cvsEpsy = self.convertNum(self.wordsNum,"#ie52s2",3,3,"float")

    def cvsEpsxy_create(self):
        self.hasAttribute("cvsName")
        self.cvsEpsxy = self.convertNum(self.wordsNum,"#ie52s2",3,4,"float")

    def cvsKapx_create(self):
        self.hasAttribute("cvsName")
        self.cvsKapx = self.convertNum(self.wordsNum,"#ie52s2",3,5,"float")

    def cvsKapy_create(self):
        self.hasAttribute("cvsName")
        self.cvsKapy = self.convertNum(self.wordsNum,"#ie52s2",3,6,"float")

    def cvsKapxy_create(self):
        self.hasAttribute("cvsName")
        self.cvsKapxy = self.convertNum(self.wordsNum,"#ie52s2",3,7,"float")

    def cvsEpsxz_create(self):
        self.hasAttribute("cvsName")
        self.cvsEpsxz = self.convertNum(self.wordsNum,"#ie52s2",3,8,"float")

    def cvsEpsyz_create(self):
        self.hasAttribute("cvsName")
        self.cvsEpsyz = self.convertNum(self.wordsNum,"#ie52s2",3,9,"float")

    def findIndice(self, att, ID):
        indice = "error"
        if "restrained" in att:
            att = "restrainedName"
        if "cbp" in att:
            att = "cbpName"
        if "jel" in att:
            att = "jelName"
        if "lnk" in att:
            att = "lnkName"
        if "gauss" in att:
            att = "gaussName"
        #print('att = {0}'.format(att))
        for i, j in enumerate(getattr(self, att)):
            if j == ID:
                indice = i
                break
        if indice == "error":
            print('indice {0} in attribute {1} does not exist'.format(ID, att))
        return indice
