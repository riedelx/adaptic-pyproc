import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# ADAPTIC class
# =============================================================================
class adaptic:

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

    def pseudoStatic(self, LF, dispMax):
        pseudoCrv = [0]
        for i in range(1,LF.shape[1]):
            pseudoTemp=np.trapz(LF[0,0:i+1],dispMax[0,0:i+1])/dispMax[0,i]
            pseudoCrv.append(pseudoTemp)
        return np.array(pseudoCrv)[None, :]

    def __init__(self,wordsNum,name):
        self.name = name
        self.wordsNum = wordsNum
        self.step = self.convertNum(wordsNum,"#io1",2,0,"int")
        self.nodeName = self.convertNum(self.wordsNum,"#in2",3,0,"str")

    def hasAttribute(self, att):
        if hasattr(self, att):
            pass
        else:
            getattr( self, str(att+"_create"))()

    def nodeDispY_create(self):
        self.nodeDispY = self.convertNum(self.wordsNum,"#in2",3,2,"float")

    def nodeDispMax_create(self):
        self.hasAttribute('nodeDispY')
        self.nodeDispMax = np.array([min(i) for i in self.nodeDispY.T])[None, :]

    def LF_create(self):
        self.LF = self.convertNum(self.wordsNum,"#io1",2,2,"float")

    def pseudo_create(self):
        self.hasAttribute('nodeDispMax')
        self.pseudo = self.pseudoStatic(self.LF, self.nodeDispMax)

    def time_create(self):
        self.time = self.convertNum(self.wordsNum,"#io1",2,2,"float")

    def gaussName_create(self):
        temp1 = self.convertNum(self.wordsNum,"#ie1s",3,0,"str")
        temp2 = self.convertNum(self.wordsNum,"#ie1s",3,1,"str")
        self.gaussName = np.array([temp1[i]+"_"+temp2[i] for i in range(len(temp2))])
        del temp1, temp2

    def restrainedName_create(self):
        self.restrainedName = self.convertNum(self.wordsNum,"#in1",3,0,"str")

    def cbpName_create(self):
        self.cbpName = self.convertNum(self.wordsNum,"#ie1",3,0,"str")

    def jelName_create(self):
        self.jelName = self.convertNum(self.wordsNum,"#ie11",3,0,"str")

    def lnkName_create(self):
        self.lnkName = self.convertNum(self.wordsNum,"#ie18",3,0,"str")

    def nodeDispX_create(self):
        self.hasAttribute("nodeName")
        self.nodeDispX = self.convertNum(self.wordsNum,"#in2",3,1,"float")

    def nodeDispY_create(self):
        self.hasAttribute("nodeName")
        self.nodeDispY = self.convertNum(self.wordsNum,"#in2",3,2,"float")

    def nodeDispRZ_create(self):
        self.hasAttribute("nodeName")
        self.nodeDispRZ = self.convertNum(self.wordsNum,"#in2",3,3,"float")

    def cbpM1_create(self):
        self.hasAttribute("cbpName")
        self.cbpM1 = self.convertNum(self.wordsNum,"#ie1",3,1,"float")

    def cbpTheta1_create(self):
        self.hasAttribute("cbpName")
        self.cbpTheta1 = self.convertNum(self.wordsNum,"#ie1d1",3,1,"float")

    def restrainedX_create(self):
        self.hasAttribute("restrainedName")
        self.restrainedX = self.convertNum(self.wordsNum,"#in1",3,1,"float")

    def restrainedY_create(self):
        self.hasAttribute("restrainedName")
        self.restrainedY = self.convertNum(self.wordsNum,"#in1",3,2,"float")

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

    def All(self):
        self.hasAttribute("nodeDispX")
        self.hasAttribute("nodeDispY")
        self.hasAttribute("nodeDispRZ")
        self.hasAttribute("cbpM1")
        self.hasAttribute("cbpM2")
        self.hasAttribute("cbpTheta1")
        self.hasAttribute("cbpTheta2")
        self.hasAttribute("cbpF")
        self.hasAttribute("cbpDelta")
        self.hasAttribute("jelF")
        self.hasAttribute("jelV")
        self.hasAttribute("jelM")
        self.hasAttribute("lnkM1")
        self.hasAttribute("lnkM2")
        self.hasAttribute("lnkF")
        self.hasAttribute("gauss1StrainB")
        self.hasAttribute("gauss1StrainT")
        self.hasAttribute("gauss2StressB")
        self.hasAttribute("gauss2StressT")
        self.hasAttribute("gauss1StrainAv")
        self.hasAttribute("gauss2StrainAv")
        self.hasAttribute("gauss1StressAv")
        self.hasAttribute("gauss2StressAv")
        self.hasAttribute("bndM1")
        self.hasAttribute("bndM2")
        self.hasAttribute("bndF")
        self.hasAttribute("bndStrain")
        self.hasAttribute("bndStress")
        self.hasAttribute("bndSlip")
        self.hasAttribute("bndBond")
        self.hasAttribute("cncM1")
        self.hasAttribute("cncM2")
        self.hasAttribute("cncF")
        self.hasAttribute("cncStrain")
        self.hasAttribute("cncStress")
        self.hasAttribute("cncGamma")
        self.hasAttribute("cncTau")

    def findIndice(self, att, ID):
        indice = "error"
        if att == "restrainedX" or att == "restrainedY" or att == "restrainedRZ":
            att = "restrainedName"
        if att == "cbpM1" or att == "cbpM2" or att == "cbpF" or att == "cbpTheta1" or att == "cbpTheta2" or att == "cbpDelta":
            att = "cbpName"
        if att == "jelF" or att == "jelV" or att == "jelM":
            att = "jelName"
        if att == "lnkM1" or att == "lnkM2" or att == "lnkF_create":
            att = "lnkName"
        if att == "gauss1StrainAv" or att == "gauss2StrainAv" or "gauss1StressAv" or att == "gauss2StressAv" or "gauss1StrainB" or att == "gauss1StrainT" or att == "gauss1StressB" or att == "gauss1StressT" or "gauss2StrainB" or att == "gauss2StrainT" or att == "gauss2StressB" or att == "gauss2StressT":
            att = "gaussName"
        if att == "nodeDispX" or att == "nodeDispY" or att == "nodeDispRZ":
            att = "nodeName"
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

    def returnVec(self,attribute,idx):
        self.hasAttribute(attribute)
        return eval('self.'+attribute+'['+str(idx)+'].T')

    # alternative constructor
    @classmethod
    def readFile(cls, title, cutoff = None,folderPath="",numPath=str('/num/')):
        wordsNum=[]
        if folderPath != "":
            path=folderPath
        else:
            path='data/'
        for data in open(path+numPath+title+str(".num"),'r'):
            wordsNum.append(data.split())
        if cutoff:
            for i in range(len(wordsNum)):
                try:
                    if wordsNum[i][0] == '#io1' and wordsNum[i+2][0] == str(cutoff):
                        break
                except:
                    pass
            wordsNum = wordsNum[0:i-2]
        return cls(wordsNum,title)
