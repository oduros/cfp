from cfp_functions import *
from numpy.linalg import eigh

# defines ground state of Er3+ ion
n = 11; l = 3; nbody = 1;
twoSplusOne1, L1, J1 = 4, 'I', 15/2
twoSplusOne2, L2, J2 = 4, 'I', 15/2

# defines the crystal field parameters
# here for Er3+ in tetragonal symmetry, C3 // z, C2 // y
# from Gaudet et al. 2017 (10.1103/PhysRevB.97.024415)
k_cef = [2,4,6]; q_cef = [-6,-3,0,3,6]
B20 = 0.0541
B40 = 0.3508
B43 = 0.1248
B60 = 0.1100
B63 = -0.0737
B66 = 0.1251

B = {2 : {0 : B20},
     4 : {-3: -B43, 0 : B40, 3 : B43},
     6 : {-6: B66, -3: -B63, 0 : B60, 3 : B63, 6 : B66}}


Jz_list = np.arange(-J1,J1+1)

hamiltonian = np.zeros((len(Jz_list),len(Jz_list)))

for k in k_cef:
    for q in q_cef:
        for iJz1, Jz1 in enumerate(np.arange(-J1, J1+1)):
            for iJz2, Jz2 in enumerate(np.arange(-J2, J2+1)):
                value = Ukq(k,q,n,l,nbody,twoSplusOne1,L1,J1,Jz1,twoSplusOne2,L2,J2,Jz2).evalf()*C_ME(k,l).evalf()
                try:
                    matrix_element = value*B[k][q]
                except:
                    matrix_element = 0
                hamiltonian[iJz1,iJz2] += matrix_element

# if one wants to print the eigenvalues (in meV)
# relative to the lowest one
eigenvalues, eigenvectors = eigh(hamiltonian)
eigenvalues = (eigenvalues - eigenvalues[0])*1000
print('Eigenvalues in meV:  ')
for eigenvalue in eigenvalues:
    print(f'{eigenvalue:03.2f}')

