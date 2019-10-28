import numpy as np
import matplotlib.pyplot as plt 
from quantum_functions import *

def wl(a):
    state_vector = matrix([[a/sqrt(1 + a^2)],[0],[0],[1/sqrt(1 + a^2)]])
    gen_state = den_mat(state_vector)
    p = var('p')
    return p*gen_state + (1 - p)*I4

def werner_touches_cvenn_at(a):
    return conditional_entropy(wl(a)).find_root(0,1)


a = np.linspace(0.01, 1, 100)

y = []
for i in a:
    y.append(werner_touches_cvenn_at(i))

plt.plot(a, y)

tnrfont = {'fontname':'Comic Sans MS'}

# plt.title("Where Werner lines touch CVENN", **tnrfont)
plt.xlabel("a (Entanglement parameter)", **tnrfont)
plt.ylabel("p (Classical mixing parameter)", **tnrfont)


plt.show()

# print(wl(1))
# print(werner_touches_cvenn_at(0.1))