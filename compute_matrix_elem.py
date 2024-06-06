from cfp_functions import *
from numpy.linalg import eig


n = 11; l = 3; nbody = 1;
twoSplusOne1, L1, J1 = 4, 'I', 15/2
twoSplusOne2, L2, J2 = 4, 'I', 15/2

A20 = 0.057024
A40 = -0.210584
A43 = 0.019893
A60 = -0.130734
A63 = 0.110532
A66 = -0.038125

A = {2 : {0 : A20},
     4 : {-3: -A43, 0 : A40, 3 : A43},
     6 : {-6: A66, -3: -A63, 0 : A60, 3 : A63, 6 : A66}}


J_list = [15/2]

hamiltonian = np.zeros((len(J_list)*len(np.arange(-max(J_list),max(J_list)+1)),(len(J_list)*len(np.arange(-max(J_list),max(J_list)+1)))))

for k in [2,4,6]:
    for q in [-6,-3,0,3,6]:
        i = 0
        for iJ1, J1 in enumerate(J_list):
            i+=len(np.arange(-J1, J1+1))
            for iJ2,J2 in enumerate(J_list):
                for iJz1, Jz1 in enumerate(np.arange(-J1, J1+1)):
                    for iJz2, Jz2 in enumerate(np.arange(-J2, J2+1)):
                        value = Ukq(k,q,n,l,nbody,twoSplusOne1,L1,J1,Jz1,twoSplusOne2,L2,J2,Jz2).evalf()*C_ME(k,l).evalf()
                        try:
                            matrix_element = value*A[k][q]
                        except:
                            matrix_element = 0
                        hamiltonian[iJ1*len(np.arange(-J1, J1+1))+iJz1,iJ2*len(np.arange(-J2, J2+1))+iJz2] += matrix_element
                
eigenvalues, eigenvectors = eig(hamiltonian)
print('Eigenvalues in meV:  ')
for eigenvalue in sorted(eigenvalues):
    egvl = eigenvalue - min(eigenvalues)
    print(f'{egvl*1000:03.15f}')