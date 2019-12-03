function x = closest_cvenn_state_2_4(rho)

    n = 8;
    
    cvx_begin 
        variable x(8,8) hermitian semidefinite
        %minimize norm(x-rho, 'fro')
        minimize trace(sqrt((x - rho)'*(x - rho)))
        %maximize trace(x*rho)
        subject to 
            trace(x) == 1
            quantum_cond_entr(x, [2 4]) >= 0
    cvx_end
    
    
    
    