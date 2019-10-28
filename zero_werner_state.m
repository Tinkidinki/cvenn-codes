function x = zero_werner_state()
    
phi_plus = [0.5 0 0 0.5; 0 0 0 0; 0 0 0 0; 0.5 0 0 0.5] 
I4 = [0.25 0 0 0; 0 0.25 0 0; 0 0 0.25 0; 0 0 0 0.25]

    n = 4;
    
    cvx_begin
        variable x(4,4) hermitian semidefinite
        variable a
        minimize quantum_cond_entr(x, [2 2]) 
        subject to 
            trace(x) == 1
            quantum_cond_entr(x, [2 2]) >= 0
            x == a*phi_plus + (1 - a)*I4
            a >= 0
            a <= 1

    cvx_end
    
    