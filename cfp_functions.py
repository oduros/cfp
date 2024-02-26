import h5py
import numpy as np
from sympy.physics.wigner import wigner_6j, wigner_3j
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy import prime, sqrt, sympify, Rational

path = './'            # path to the CFP database
cfp_db=h5py.File(f'{path}cfp.h5','r')

def readLSalpha(term): #returns L, S and if applicable, J for a given spectroscopic term
    if len(term)==2:
        S,Lstr = int(term[0]),term[1]
        alpha = '0'
    else:
        S,Lstr,alpha = int(term[0]),term[1],float(term[2:])
    if   Lstr == 'S': L=0
    elif Lstr == 'P': L=1
    elif Lstr == 'D': L=2
    elif Lstr == 'F': L=3
    elif Lstr == 'G': L=4
    elif Lstr == 'H': L=5
    elif Lstr == 'I': L=6
    elif Lstr == 'J': L=7
    elif Lstr == 'K': L=8
    elif Lstr == 'L': L=9
    elif Lstr == 'M': L=10
    elif Lstr == 'N': L=11
    elif Lstr == 'O': L=12
    elif Lstr == 'P': L=13
    elif Lstr == 'Q': L=14
    elif Lstr == 'R': L=15
    elif Lstr == 'S': L=16
    elif Lstr == 'T': L=17
    elif Lstr == 'U': L=18
    elif Lstr == 'V': L=19
    elif Lstr == 'W': L=20
    elif Lstr == 'X': L=21
    elif Lstr == 'Y': L=22
    elif Lstr == 'Z': L=23
    else: print('No L value found for term',term)
    return L,(S-1)/2, alpha
def readpair(pair):
    coma = pair.index(',')
    parent1 = pair[0:coma]
    parent2 = pair[coma+1:]
    return parent1, parent2

def read_seniority(conf, term):
    return cfp_db['states'][f'{conf}'][f'{term}']['seniority_v'][()]

def phase_bf(v):
    return 1 if v in [0,1,4,5,8,9,12,13,16,17] else -1

def U_matrix_element(k, n, l, nbody, twoSplusOne1, L1, twoSplusOne2, L2, alpha1, alpha2):
    term, term_prime = f'{twoSplusOne1}{L1}{alpha1}', f'{twoSplusOne2}{L2}{alpha2}' #creates strings to search for the two desired terms of the l^n configuration in the CFP database
    daughter_term, daughter_term_prime = term, term_prime
    # c1, c2 = read_seniority(conf, term), read_seniority(conf, term_prime)
    L, S, alpha  = readLSalpha(daughter_term)                            #converts the input parameters of the spectroscopic term ^(2S+1)L_alpha into the associated quantum numbers
    Lprime, Sprime, alphaprime = readLSalpha(daughter_term_prime)       #same for the second state
    if n > 2 * l + 1 :                                              #checks whether the l orbital is more than half-filled or not
        w = l * 4 + 2 - n                                           #if yes, follows Cowan 11.58 definition to compute the matrix element of the term from the less-than-half-filled configuration equivalent
        more_than_half_filled = True
    else:
        w = n
        more_than_half_filled = False
    if l == 0:                                                      #creates string to search for the correct configuration in the CFP database
        conf = f's{w}'
    elif l == 1:
        conf = f'p{w}'
    elif l == 2:
        conf = f'd{w}'
    elif l == 3:
        conf = f'f{w}'
    else:
        print(f'CFP for {l}-type electrons are not yet listed in the dataset.');

    if w == 1:
        matrix_element = int(1)                                     #when filling is w = 1 (1 electron or 1 hole in the orbital), the matrix element is reduced to 1
         
    else:
        sum = sympify(0)                                            #creating the sum on which will be summed the computations of CFP, following Cowan 11.53
                                                                    #Using sympify() from SymPy to be sure the output will be a symbolic quantity
        for parents in cfp_db[f'{nbody}'][f'{conf}'][daughter_term].keys():  # reading over all CFPs of the daughter_term
            parent1, parent2 = readpair(parents)
            Lsum, Ssum, alphasum = readLSalpha(parent1)                     #defining Lsum, Ssum and alphasum, the quantum numbers used for the summing  in Cowan 11.53
            CFPstring1 = cfp_db[f'{nbody}'][f'{conf}'][daughter_term][parents][:] #Extracting the CFP as the string it is stored as in the database.
            CFP1 = sympify(int(CFPstring1[0]))                      #The string has always at least one parameter (the a0 definded by Nielson & Koster)
            if len(CFPstring1) > 1:                                 #Then, if the CFP is longer than only a0, we will read a1, a2...an where n is the lenght of the CFPstring
                for n in range(len(CFPstring1) - 1):
                    CFP1 *= prime(n + 1) ** sympify(int(CFPstring1[n + 1]))     #We multiply a0 by the first prime of the prime series [2,3,5,7,...,29,31] to the power of a1
            else:
                CFP1 = CFP1                                         #If CFPstring is no longer than 1, our CFP is just a0

            sign = np.sign(CFP1)                                    #Reading and storing the sign of our calculated CFP
            CFP1 = sign * sqrt(np.abs(CFP1))                        #Normalization of Cowan 9.43 requires that the final CFP is its sign times its square root (see also Cowan 9.45)

            try:                                                    #We use a Try and Except statement to check if, for the second term of the configuration called, a CFP exists with the same parents
                CFPstring2 = cfp_db[f'{nbody}'][f'{conf}'][daughter_term_prime][parents][:]       #Then, same calculation is done for CFP2 as it was for CFP1
                CFP2 = sympify(int(CFPstring2[0]))
                if len(CFPstring2) > 1:
                    for n in range(len(CFPstring2) - 1):
                        CFP2 *= prime(n + 1) ** sympify(int(CFPstring2[n + 1]))
                else:
                    CFP2 = CFP2
                sign = np.sign(CFP2)
                CFP2 = sign * sqrt(np.abs(CFP2))
            except:                                                 #If no CFP with the same parents exists for the second term, we just continue without any CFP2
                continue
            try:                                                    #Try and Except statement as we might not have any CFP2, or a null value of the 6j symbol
                sum += sympify((-1) ** (Lsum)) * wigner_6j(l, k, l, L, Lsum, Lprime) * (CFP1 * CFP2)  #Summing following Cowan 11.53. wigner_6j() returns a symbolic quantity
            except:
                continue
        matrix_element = KroneckerDelta(S, Sprime) * sympify(w) * sympify((-1) ** (l + L + k)) * sqrt(
            (int(2 * L) + 1) * (int(2 * Lprime) + 1),evaluate = False) * sum        #Cowan 11.53, "part two", once the sum has been calculated. The parameter 'evaluate' takes a bool and allows
                                                                                    #to be sure that we stay in symbolic quantities and not floats.
    if more_than_half_filled:
        matrix_element = matrix_element * (1 / sympify(-(-1) ** k))             #Cowan 11.58 for more than half-filled shells
    else:                                                                       #If not more than half-filled, the value does not change
        matrix_element = matrix_element
    return matrix_element, L, S, Lprime, Sprime                                     #The matrix element is here returned along L, S, L' and S' for later purpose, especially in function Ukq() defined below   


def Ukq(k, q, n, l, nbody, twoSplusOne1, L1, J1, Jz1, twoSplusOne2, L2, J2, Jz2, alpha1 = '', alpha2 = ''):         #alpha1 and alpha2 are optional parameters for some configurations but mandatory 
                                                                                                                    #to determine the cfp to use when more than one term is possible for a pair of L and S
    matrix_element, L, S, Lprime, Sprime = U_matrix_element(k, n, l, nbody, twoSplusOne1, L1, twoSplusOne2, L2,alpha1,alpha2)    #loading the matrix element for the desired terms and configuration
    J = J1 ; Jprime = J2 ; Jz = Jz1 ; Jzprime = Jz2                #homogeneous writing of the quantum numbers
    Uk = ((-1)**Rational(S+L+Jprime+k))*sqrt((int(2*J)+1)*(int(2*Jprime)+1),evaluate = False)*wigner_6j(J, Jprime,k, Lprime,L,S)* matrix_element #Adding J, following Wybourne 6-5
    Ukq = ((-1)**Rational(J-Jz))*wigner_3j(J,k,Jprime,-Jz,q,Jzprime)*Uk         #Adding Jz, following Wybourne 6-4
    return Ukq

def Ukq_float(k, q, n, l, nbody, twoSplusOne1, L1, J1, Jz1, twoSplusOne2, L2, J2, Jz2, alpha1 = '', alpha2 = ''):
    Ukq_ft = Ukq(k, q, n, l, nbody, twoSplusOne1, L1, J1, Jz1, twoSplusOne2, L2, J2, Jz2, alpha1, alpha2)
    return Ukq_ft.evalf()