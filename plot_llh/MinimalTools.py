import numpy as np
from PhysConst import PhysicsConstants

def eigenvectors(M):
        """ Calculates the eigenvectors and eigenvalues ordered by eigenvalue size

        @type  M     : matrix
        @param M     : matrix M

        @rtype       : list
        @return      : [eigenvalues list, eigenvector list]
        """
        D,V = np.linalg.eig(M)
        DV = []
        VT = V.T
        for i,eigenvalue in enumerate(D):
            DV.append([eigenvalue,VT[i]]) 

        DV = sorted(DV,key = lambda x : x[0].real)#np.abs(x[0].real))

        V2 = []
        D2 = []
        for e in DV:
            V2.append(e[1])
            D2.append(e[0])
        return D2,V2

# General Rotation Matrix
def R(i,j,cp,param):
        """ Rotation Matrix
        Calculates the R_ij rotations. Also incorporates CP-phases when necesary.
        @type   i       :   int
        @param  i       :   i-column.
        @type   j       :   int
        @param  j       :   j-row.
        @type   cp      :   int
        @param  cp      :   if cp = 0 : no CP-phase. else CP-phase = CP_array[cp]

        @rtype          :   numpy.array
        @return         :   returns the R_ij rotation matrix.
        """
        # if cp = 0 -> no complex phase
        # R_ij, i<j
        if(j<i):
            # no funny business
            l = i
            i = j
            j = l

        # rotation matrix generator
        R = np.zeros([param.numneu,param.numneu],complex)
        # diagonal terms
        for k in np.arange(0,param.numneu,1):
            if(k != i-1 and k != j-1):
                R[k,k] = 1.0
            else :
                R[k,k] = param.c[i,j]
        # non-diagonal terms
        if(cp != 0):
            sd = np.sin(param.dcp[cp])
            cd = np.cos(param.dcp[cp])
            faseCP = complex(cd,sd)
        else:
            faseCP = complex(1.0,0.0)

        R[i-1,j-1] = param.s[i,j]*faseCP.conjugate()
        R[j-1,i-1] = -param.s[i,j]*faseCP
        return R

def calcU(param):
        """ Defining the mixing matrix parametrization.
        @type   param   :   PhysicsConstants
        @param  param   :   set of physical parameters to be used.

        @rtype          :   None
        @return         :   Sets mixing matrix.
        """
        if(param.numneu == 2):
	        return self.R(1,2,0,param)
        elif(param.numneu == 3):
            return np.dot(R(2,3,0,param),np.dot(R(1,3,1,param),R(1,2,0,param)))
        elif(param.numneu == 4):
            return np.dot(R(3,4,0,param),np.dot(R(2,4,2,param),np.dot(R(1,4,0,param),np.dot(R(2,3,0,param),np.dot(R(1,3,1,param),R(1,2,0,param))))))
        elif(param.numneu == 5):
            return np.dot(R(4,5,0,param),np.dot(R(3,5,0,param),np.dot(R(2,5,0,param),np.dot(R(1,5,3,param),np.dot(R(3,4,0,param),np.dot(R(2,4,0,param),np.dot(R(1,4,2,param),np.dot(R(2,3,0,param),np.dot(R(1,3,1,param),R(1,2,0,param))))))))))
        elif(param.numneu == 6):
            # 3+3 twin-sterile-neutrino model
            return np.dot(R(3,6,0,param),np.dot(R(2,5,0,param),np.dot(R(1,4,0,param),np.dot(R(2,3,0,param),np.dot(R(1,3,1,param),R(1,2,0,param))))))
        else:
            print "Sorry, too many neutrinos. Not yet implemented! =(."
            quit()

        # antineutrino case
        if param.neutype == "antineutrino" :
            return self.U.conjugate()



def massM2(param):
    """ Mass term in the neutrino mass basis.

    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.

    @rtype          :   numpy array
    @return         :   mass matrix in mass basis.
    """
    M2 = np.zeros([param.numneu,param.numneu],complex)
    for k in np.arange(1,param.numneu,1):
        M2[k,k] = param.dmsq[k+1]
    return M2

def flavorM2(param):
    """ Mass term in the neutrino flavor basis.

    @type   param   :   PhysicsConstants
    @param  param   :   set of physical parameters to be used.

    @rtype          :   numpy array
    @return         :   mass matrix in flavor basis.
    """
    U = calcU(param)
    return np.dot(U,np.dot(massM2(param),U.conjugate().T))

class LVP(PhysicsConstants):
    def __init__(self):
        super(LVP,self).__init__()
        self.th12 = 0.0
        self.th13 = 0.0
        self.th23 = 0.0
        self.delta1 = 0.0
        self.deltaCP = 0.0
        super(LVP,self).Refresh()

        # LVS
        self.LVS21 = 0.0	                #
        self.LVS31 = 0.0                    #
        self.LVS41 = 0.0                    #
        self.LVS51 = 0.0                    #
        self.LVS61 = 0.0                    #
        # SQUARED MASS DIFFERENCE MATRIX
        self.LVS = np.zeros([self.numneumax+2],float)
        self.LVS[2] = self.LVS21
        self.LVS[3] = self.LVS31
        self.LVS[4] = self.LVS41
        self.LVS[5] = self.LVS51
        self.LVS[6] = self.LVS61

    def Refresh(self):
        super(LVP,self).Refresh()
        LVS = self.LVS
        LVS[2] = self.LVS21
        LVS[3] = self.LVS31
        LVS[4] = self.LVS41
        LVS[5] = self.LVS51
        LVS[6] = self.LVS61

def DiagonalMatrixLV(param):
    DD = np.zeros([param.numneu,param.numneu],complex)
    for k in np.arange(1,param.numneu,1):
        DD[k,k] = param.LVS[k+1]
    return DD

def LVTerm(LVparam):
    U = calcU(LVparam)
    return np.dot(U,np.dot(DiagonalMatrixLV(LVparam),U.conjugate().T))
