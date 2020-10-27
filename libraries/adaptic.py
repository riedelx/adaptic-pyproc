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
                        tempVal = float(wordsNum[jRow][column])
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

    def gaussName_create(self):
        temp1 = self.convertNum(self.wordsNum,"#ie1s",3,0,"str")
        temp2 = self.convertNum(self.wordsNum,"#ie1s",3,1,"str")
        self.gaussName = np.array([temp1[i]+"_"+temp2[i] for i in range(len(temp2))])
        del temp1, temp2

    def restrainedName_create(self):
        self.restrainedName = self.convertNum(self.wordsNum,"#in1",3,0,"str")

    def elementName_create(self):
        self.elementName = self.convertNum(self.wordsNum,"#ie1",3,0,"str")

    def jelName_create(self):
        self.jelName = self.convertNum(self.wordsNum,"#ie11",3,0,"str")

    def lnk2Name_create(self):
        self.lnk2Name = self.convertNum(self.wordsNum,"#ie18",3,0,"str")

    def nodeDispX_create(self):
        self.hasAttribute("nodeName")
        self.nodeDispX = self.convertNum(self.wordsNum,"#in2",3,1,"float")

    def nodeDispY_create(self):
        self.hasAttribute("nodeName")
        self.nodeDispY = self.convertNum(self.wordsNum,"#in2",3,2,"float")

    def nodeDispRZ_create(self):
        self.hasAttribute("nodeName")
        self.nodeDispRZ = self.convertNum(self.wordsNum,"#in2",3,3,"float")

    def elementM1_create(self):
        self.hasAttribute("elementName")
        self.elementM1 = self.convertNum(self.wordsNum,"#ie1",3,1,"float")

    def restrainedX_create(self):
        self.hasAttribute("restrainedName")
        self.restrainedX = self.convertNum(self.wordsNum,"#in1",3,1,"float")

    def restrainedY_create(self):
        self.hasAttribute("restrainedName")
        self.restrainedY = self.convertNum(self.wordsNum,"#in1",3,2,"float")

    def restrainedRZ_create(self):
        self.hasAttribute("restrainedName")
        self.restrainedRZ = self.convertNum(self.wordsNum,"#in1",3,3,"float")

    def elementM2_create(self):
        self.hasAttribute("elementName")
        self.elementM2 = self.convertNum(self.wordsNum,"#ie1",3,2,"float")

    def elementF_create(self):
        self.hasAttribute("elementName")
        self.elementF = self.convertNum(self.wordsNum,"#ie1",3,3,"float")

    def jelF_create(self):
        self.hasAttribute("jelName")
        self.jelF = self.convertNum(self.wordsNum,"#ie11",3,1,"float")

    def jelV_create(self):
        self.hasAttribute("jelName")
        self.jelV = self.convertNum(self.wordsNum,"#ie11",3,2,"float")

    def jelM_create(self):
        self.hasAttribute("jelName")
        self.jelM = self.convertNum(self.wordsNum,"#ie11",3,3,"float")

    def lnk2M1_create(self):
        self.hasAttribute("lnk2Name")
        self.lnk2M1 = self.convertNum(self.wordsNum,"#ie18",3,1,"float")

    def lnk2M2_create(self):
        self.hasAttribute("lnk2Name")
        self.lnk2M2 = self.convertNum(self.wordsNum,"#ie18",3,2,"float")

    def lnk2F_create(self):
        self.hasAttribute("lnk2Name")
        self.lnk2F = self.convertNum(self.wordsNum,"#ie18",3,3,"float")

    def gauss1StrainB_create(self):
        self.hasAttribute("gaussName")
        self.gauss1StrainB = self.convertNum(self.wordsNum,"#ie1s",3,2,"float")

    def gauss1StrainT_create(self):
        self.hasAttribute("gaussName")
        self.gauss1StrainT = self.convertNum(self.wordsNum,"#ie1s",3,4,"float")

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

    def All(self):
        self.hasAttribute("nodeDispX")
        self.hasAttribute("nodeDispY")
        self.hasAttribute("nodeDispRZ")
        self.hasAttribute("elementM1")
        self.hasAttribute("elementM2")
        self.hasAttribute("elementF")
        self.hasAttribute("jelF")
        self.hasAttribute("jelV")
        self.hasAttribute("jelM")
        self.hasAttribute("lnk2M1")
        self.hasAttribute("lnk2M2")
        self.hasAttribute("lnk2F")
        self.hasAttribute("gauss1StrainB")
        self.hasAttribute("gauss1StrainT")
        self.hasAttribute("gauss2StressB")
        self.hasAttribute("gauss2StressT")

    def findIndice(self, att, ID):
        indice = "error"
        if att == "restrainedX" or att == "restrainedY" or att == "restrainedRZ":
            att = "restrainedName"
        if att == "elementM1" or att == "elementM2" or att == "elementF":
            att = "elementName"
        if att == "jelF" or att == "jelV" or att == "jelM":
            att = "jelName"
        if att == "lnk2M1" or att == "lnk2M2" or att == "lnk2F_create":
            att = "lnk2Name"
        if att == "gauss1StrainB" or att == "gauss1StrainT" or att == "gauss1StressB" or att == "gauss1StressT" or "gauss2StrainB" or att == "gauss2StrainT" or att == "gauss2StressB" or att == "gauss2StressT":
            att = "gaussName"
        if att == "nodeDispX" or att == "nodeDispY" or att == "nodeDispRZ":
            att = "nodeName"
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
    def readFile(cls, title, cutoff = None,folderPath=""):
        wordsNum=[]
        if folderPath != "":
            path=folderPath
        else:
            path='data/'
        for data in open(path+str('/num/')+title+str(".num"),'r'):
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
