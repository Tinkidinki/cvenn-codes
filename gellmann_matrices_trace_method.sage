
lambda_1 = matrix([[0, 1, 0],[1, 0, 0],[0, 0, 0]], ring=CC)/sqrt(2)
lambda_2 = matrix([[0, -i, 0],[i, 0, 0],[0, 0, 0]], ring=CC)/sqrt(2)

lambda_3 = matrix([[1, 0, 0],[0, -1, 0],[0, 0, 0]], ring=CC)/sqrt(2)

lambda_4 = matrix([[0, 0, 1],[0, 0, 0],[1, 0, 0]], ring=CC)/sqrt(2)
lambda_5 = matrix([[0, 0, -i],[0, 0, 0],[i, 0, 0]], ring=CC)/sqrt(2)

lambda_6 = matrix([[0, 0, 0],[0, 0, 1],[0, 1, 0]], ring=CC)/sqrt(2)
lambda_7 = matrix([[0, 0, 0],[0, 0, -i],[0, i, 0]], ring=CC)/sqrt(2)

lambda_8 = (1/3^(1/2)) * matrix([[1, 0, 0],[0, 1, 0],[0, 0, -2]], ring=CC)/sqrt(2)
I_3 = matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]], ring=CC)/sqrt(3)

gb = [lambda_1, lambda_2, lambda_3, 
    lambda_4, lambda_5, lambda_6, lambda_7, lambda_8, I_3]

def convert_to_gellmann_basis(rho):
    gellmann_matrix = matrix([[0 for i in range(9)] for i in range(9)], ring=CC)
    for p in range(len(gb)):
        for q in range(len(gb)):
            gellmann_matrix[p,q] = (rho*gb[p].tensor_product(gb[q])).trace()
    return gellmann_matrix

def convert_to_comp_basis(rho):
    comp_matrix = matrix([[0 for i in range(9)] for i in range(9)], ring=CC)
    for i in range(9):
        for j in range(9):
            comp_matrix = comp_matrix + rho[i,j]*gb[i].tensor_product(gb[j])
    return comp_matrix

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

