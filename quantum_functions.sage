
from variables import VariableGenerator
import numpy

# Pauli Matrices
X = Matrix(QQbar, 2, 2, [0, 1, 1, 0])
Y = Matrix(QQbar, 2, 2, [0, -I, I, 0])
Z = Matrix(QQbar, 2, 2, [1, 0, 0, -1])

Pauli = [X,Y,Z]
I2 = Matrix(QQbar, 2, 2, [1/2, 0, 0, 1/2])
I3 = Matrix(QQbar, 3, 3, [1/3, 0, 0, 0, 1/3, 0, 0, 0, 1/3])
I4 = Matrix(QQbar, 4, 4, [1/4, 0, 0, 0, 0, 1/4, 0, 0, 0, 0, 1/4, 0, 0, 0, 0, 1/4])
Identity4 = Matrix(QQbar, 4, 4, [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
Identity3 = Matrix(QQbar, 3, 3, [1, 0, 0, 0, 1, 0, 0, 0, 1])
Identity2 = Matrix(QQbar, 2, 2, [1, 0, 0, 1])

phi_plus = matrix([[1/sqrt(2)],[0],[0],[1/sqrt(2)]])
phi_minus = matrix([[1/sqrt(2)],[0],[0],[-1/sqrt(2)]])
psi_plus = matrix([[0],[1/sqrt(2)],[1/sqrt(2)],[0]])
psi_minus = matrix([[0],[1/sqrt(2)],[-1/sqrt(2)],[0]])

def den_mat(v):
    return v*v.conjugate_transpose()

def entropy(M):
    eigs = M.eigenvalues()
    entropy = 0.0
    for i in eigs:
        if (i!=0):
            entropy += -i*log(i,2)
    return entropy

def trace(M):
    trace = 0
    for i in range(len(M.rows())):
        trace += M[i][i]
    return trace

def partial_trace(Q):
    return  Matrix(SR, 2, 2, [Q[0][0]+Q[1][1], Q[0][2]+Q[1][3], Q[2][0]+Q[3][1], Q[2][2]+Q[3][3]])

def conditional_entropy(M):
    return entropy(M) - entropy(partial_trace(M))

def Bell_diag():
    a = VariableGenerator('a')
    return Matrix(SR, 4, 4, [a[1] + a[2], 0, 0, a[1] - a[2], 0, a[3]+a[4], a[3]-a[4], 0, 0, a[3]-a[4], a[3]+a[4], 0, a[1]-a[2], 0, 0, a[1]+a[2]])/2



def bloch_r(M):
    r = [0, 0, 0]
    for i in range(3):
        r[i] = trace(M * Pauli[i].tensor_product(I2))
    return r

def bloch_s(M):
    s = [0, 0, 0]
    for i in range(3):
        s[i] = trace(M * I2.tensor_product(Pauli[i]))
    return s

def bloch_t(M):
    t = Matrix(M.base_ring(), 3, 3, [0, 0, 0, 0, 0, 0, 0, 0, 0])
    for i in range(3):
        for j in range(3):
            t[i, j] = trace(M * Pauli[i].tensor_product(Pauli[j]))
    return t

def square_root(M):
    D, S = M.eigenmatrix_right()
    return S*diagonal_matrix([sqrt(x) for x in D.diagonal()]) *S^-1

def partial_transpose(M):
# Need a better way of writing this, and need to make the method more general.

    array = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]

    array[0][0] = M[0,0]
    array[0][1] = M[1,0]
    array[1][0] = M[0,1]
    array[1][1] = M[1,1]
    array[0][2] = M[0,2]
    array[0][3] = M[1,2]
    array[1][2] = M[0,3]
    array[1][3] = M[1,3]
    array[2][0] = M[2,0]
    array[2][1] = M[3,0]
    array[3][0] = M[2,1]
    array[3][1] = M[3,1]
    array[2][2] = M[2,2]
    array[2][3] = M[3,2]
    array[3][2] = M[2,3]
    array[3][3] = M[3,3]

    return Matrix(array)

def bloch_to_matrix(bloch_r, bloch_s, bloch_T):
    M = matrix(4, ring=QQ)
   
    for i in range(3):
        M += bloch_r[i,0]*(Pauli[i].tensor_product(Identity2))
        print M
        M += bloch_s[i,0]*(Identity2.tensor_product(Pauli[i]))
        print M

    for i in range(3):
        for j in range(3):
            M += bloch_T[i,j]*(Pauli[i].tensor_product(Pauli[j]))
            print M

    M += Identity4
    M/=4

    return M

def matrix_to_bloch(M):
    return (bloch_r(M), bloch_s(M), bloch_t(M))

def peres_horodecki(M):

    w3_mat = matrix([[M[0,0], M[1,0],M[0,2]],
    [M[0,1],M[1,1],M[0,3]],
    [M[2,0], M[3,0], M[2,2]]])

    w4_mat = matrix([[M[0,0],M[1,0],M[0,2],M[1,2]],
    [M[0,1],M[1,1],M[0,3],M[1,3]],
    [M[2,0],M[3,0],M[2,2],M[3,2]],
    [M[2,1], M[3,1], M[2,3], M[3,3]]])

    w2_mat = matrix([[M[0,0],M[1,0]],[M[0,1],M[1,1]]])

    w3 = w3_mat.determinant()
    w4 = w4_mat.determinant()
    w2 = w2_mat.determinant()

    print "w3_mat\n", w3_mat
    print "w4_mat\n", w4_mat
    print "w2_mat\n", w2_mat
    print w3
    print w4
    print w2

    a = var('a')
    b = var('b')
    p = var('p')

    # print "w3<0"
    # print solve(w3<0, a, b, p)
    # print "w4<0"
    # print solve(w3<0, a, b, p)
    # print "w2>=0"
    # print solve(w2>=0, a, b, p)


    # if ((w2>=0) and (w3<0 or w4<0)):
    #     return 1
    # else:
    #     return 0

def test():
    print "test"

def rho():
    a = VariableGenerator('a')
    m = matrix([[a[0], a[4] + I*a[5], a[6] + I*a[7], a[8] + I*a[9]],[a[4] - I*a[5], a[1], a[10] + I*a[11], a[12]+I*a[13]],[a[6] - I*a[7], a[10] - I*a[11], a[2], a[14] + I*a[15]],[a[8] - I*a[9], a[12] - I*a[13], a[14] - I*a[15], a[3]]])
    return m




