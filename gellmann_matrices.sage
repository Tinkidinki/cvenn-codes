
#Convention followed for order:
#Symmetric matrices
#Anti symmetric matrices
#Diagonal matrices
from sage.misc.decorators import infix_operator




lambda_1 = matrix([[0, 1, 0],[1, 0, 0],[0, 0, 0]])
lambda_2 = matrix([[0, -i, 0],[i, 0, 0],[0, 0, 0]])

lambda_3 = matrix([[1, 0, 0],[0, -1, 0],[0, 0, 0]])

lambda_4 = matrix([[0, 0, 1],[0, 0, 0],[1, 0, 0]])
lambda_5 = matrix([[0, 0, -i],[0, 0, 0],[i, 0, 0]])

lambda_6 = matrix([[0, 0, 0],[0, 0, 1],[0, 1, 0]])
lambda_7 = matrix([[0, 0, 0],[0, 0, -i],[0, i, 0]])

lambda_8 = (1/3^(1/2)) * matrix([[1, 0, 0],[0, 1, 0],[0, 0, -2]])
I = matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])

b00 = var('b00')
b01 = var('b01')
b02 = var('b02')

b10 = var('b10')
b11 = var('b11')
b12 = var('b12')

b20 = var('b20')
b21 = var('b21')
b22 = var('b22')

comp_basis = [b00, b01, b02, b10, b11, b12, b20, b21, b22]

l1 = var('l1')
l2 = var('l2')
l3 = var('l3')
l4 = var('l4')
l5 = var('l5')
l6 = var('l6')
l7 = var('l7')
l8 = var('l8')
I_3 = var('I_3')

gellmann_basis = [l1, l2, l3, l4, l5, l6, l7, l8, I_3]

# R.<b00, b01, b02, b10, b11, b12, b20, b21, b22, l1, l2, l3, l4, l5, l6, l7, l8, I_3, i> = FreeAlgebra(CC)

def gellmann_eq(v):
    return {
        var('b00'):(3*l3 + sqrt(3).n()*l8 + 2*I_3)/6,
        var('b01'):(l1 + i*l2)/2,
        var('b02'):(l4 + i*l5)/2,
        var('b10'):(l1 - i*l2)/2,
        var('b11'):(-3*l3 + sqrt(3).n()*l8 + 2*I_3)/6,
        var('b12'):(l6 + i*l7)/2,
        var('b20'):(l4 - i*l5)/2,
        var('b21'):(l6 - i*l7)/2,
        var('b22'):(I_3 - sqrt(3).n()*l8)/3
    }[v]

def comp_basis_eq(v):
    return {
        var('l1') : b01 + b10,
        var('l2') : -i*b01 + i*b10,
        var('l3') : b00 - b11,
        var('l4') : b02 + b20,
        var('l5') : -i*b02 + i*b20,
        var('l6') : b12 + b21,
        var('l7') : -i*b12 + i*b21,
        var('l8') : (1/3^(1/2))*(b00 + b11 -2*b22),
        var('I_3'): b00 + b11 + b22,
    }[v]

def convert_to_gellmann(expr):
   substituted_expr = expr.subs({var : gellmann_eq(var) for var in comp_basis})
   return substituted_expr.full_simplify()

def convert_to_comp_basis(expr):
    substituted_expr = expr.subs({var: comp_basis_eq(var) for var in gellmann_basis})
    return substituted_expr.full_simplify()

##### Testing purposes: analytical witness ####
aw = 0.0721*b00*b00 + 0.1306*b01*b01 + 0.1306*b02*b02 + 0.1306*b10*b10 + 0.0721*b11*b11 + 0.1306*b12*b12 + 0.1306*b20*b20 + 0.1306*b21*b21 + 0.0721*b22*b22 -0.0584 *(b11*b00 + b22*b00 + b00*b11 + b00*b22 + b11*b22 + b22*b11)
