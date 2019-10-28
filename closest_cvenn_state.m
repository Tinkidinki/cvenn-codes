function x = closest_cvenn_state(rho)

n = 4;

cvx_begin 
    variable x(4,4) hermitian semidefinite
    minimize norm(x-rho, 'fro')
    %maximize trace(x*rho)
    subject to 
        trace(x) == 1
        quantum_cond_entr(x, [2 2]) >= 0 
cvx_end



