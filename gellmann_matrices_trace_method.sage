'''
Change of basis to a non-normalised orthogonal basis:
We have two options when we do this.
Option 1: Change to normalised basis, and then divide by normalisation constant.
So, a) Normalise the basis
    b) Take trace (or inner product) of object in computational basis with new basis to find coefficients
    c) Divide coefficients by the normalisation constant. 

Option 2: We know in the process that there are two extra multiplications by normalisation constant. So, 
    a) Take inner product with respect to new basis only
    b) Divide by square of normalisation constant to compensate for the double multiplication. 

Option 2 is being done here
'''

''' The Computational Basis'''
c00 = matrix([[1, 0, 0],[0, 0, 0],[0, 0, 0]])
c01 = matrix([[0, 1, 0],[0, 0, 0],[0, 0, 0]])
c02 = matrix([[0, 0, 1],[0, 0, 0],[0, 0, 0]])

c10 = matrix([[0, 0, 0],[1, 0, 0],[0, 0, 0]])
c11 = matrix([[0, 0, 0],[0, 1, 0],[0, 0, 0]])
c12 = matrix([[0, 0, 0],[0, 0, 1],[0, 0, 0]])

c20 = matrix([[0, 0, 0],[0, 0, 0],[1, 0, 0]])
c21 = matrix([[0, 0, 0],[0, 0, 0],[0, 1, 0]])
c22 = matrix([[0, 0, 0],[0, 0, 0],[0, 0, 1]])


'''The Gellmann Matrices'''
lambda_1 = matrix([[0, 1, 0],[1, 0, 0],[0, 0, 0]], ring=CC)
lambda_2 = matrix([[0, -i, 0],[i, 0, 0],[0, 0, 0]], ring=CC)

lambda_3 = matrix([[1, 0, 0],[0, -1, 0],[0, 0, 0]], ring=CC)

lambda_4 = matrix([[0, 0, 1],[0, 0, 0],[1, 0, 0]], ring=CC)
lambda_5 = matrix([[0, 0, -i],[0, 0, 0],[i, 0, 0]], ring=CC)

lambda_6 = matrix([[0, 0, 0],[0, 0, 1],[0, 1, 0]], ring=CC)
lambda_7 = matrix([[0, 0, 0],[0, 0, -i],[0, i, 0]], ring=CC)

lambda_8 = (1/3^(1/2)) * matrix([[1, 0, 0],[0, 1, 0],[0, 0, -2]], ring=CC)
I_3 = matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]], ring=CC)

'''Gellmann Basis'''
gb = [lambda_1, lambda_2, lambda_3, 
    lambda_4, lambda_5, lambda_6, lambda_7, lambda_8, I_3]

'''Gellmann Basis divided by square of normalisation constant'''
gb_n = [lambda_1/2, lambda_2/2, lambda_3/2, 
    lambda_4/2, lambda_5/2, lambda_6/2, lambda_7/2, lambda_8/2, I_3/3]

def convert_to_gellmann_basis(rho):
    gellmann_matrix = matrix([[0 for i in range(9)] for i in range(9)], ring=CC)
    for p in range(len(gb)):
        for q in range(len(gb)):
            gellmann_matrix[p,q] = (rho*gb_n[p].tensor_product(gb_n[q])).trace()
    return gellmann_matrix.n(20)

def convert_to_comp_basis(rho):
    comp_matrix = matrix([[0 for i in range(9)] for i in range(9)], ring=CC)
    for i in range(9):
        for j in range(9):
            comp_matrix = comp_matrix + rho[i,j]*gb[i].tensor_product(gb[j])
    return comp_matrix.n(20)

''' Test matrix'''
w = matrix([[ 0.0721,         0,         0,         0,   -0.0584,    0.0000,    0.0000,   -0.0000,   -0.0584],
            [      0,    0.1306,         0,         0,         0,         0,         0,         0,         0],
            [      0,         0,    0.1306,         0,         0,         0,         0,         0,         0],
            [      0,         0,         0,    0.1306,         0,         0,         0,         0,         0],
            [-0.0584,         0,         0,         0,    0.0721,   -0.0000,   -0.0000,    0.0000,   -0.0584],
            [ 0.0000,         0,         0,         0,   -0.0000,    0.1306,    0.0000,    0.0000,   -0.0000],
            [ 0.0000,         0,         0,         0,   -0.0000,    0.0000,    0.1306,    0.0000,   -0.0000],
            [-0.0000,         0,         0,         0,    0.0000,    0.0000,    0.0000,    0.1306,    0.0000],
            [-0.0584,         0,         0,         0,   -0.0584,   -0.0000,   -0.0000,    0.0000,    0.0721]]\
            , ring=CC)

def get_witness(a, b, c):
    return a* (c00.tensor_product(c00) + c11.tensor_product(c11) + c22.tensor_product(c22))+ \
    b*(c00.tensor_product(c11) + c00.tensor_product(c22) + c11.tensor_product(c00) + c11.tensor_product(c22) + c22.tensor_product(c00) + c22.tensor_product(c11)) + \
    c*(c01.tensor_product(c01) + c02.tensor_product(c02) + c10.tensor_product(c10) + c12.tensor_product(c12) + c20.tensor_product(c20) + c21.tensor_product(c21))
    