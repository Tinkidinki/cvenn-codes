function x = closest_cvenn_state_3D(rho)

    n = 9;
    
    cvx_begin 
        variable x(9,9) hermitian semidefinite
        minimize norm(x-rho, 'fro')
        %maximize trace(x*rho)
        subject to 
            trace(x) == 1
            quantum_cond_entr(x, [3 3]) >= 0
    cvx_end
    
    
    
    