function cvx_optval = tight_upper_bound(sigma, chi)

n = 4;

cvx_begin
    variable rho(4,4) hermitian semidefinite
    maximize quantum_cond_entr(rho, [2 2]) / log(2)
    
    subject to 
        trace(rho) == 1

        -trace(rho * logm(sigma)/log(2)) + trace(TrX(rho, [1], [2,2]) * logm(TrX(sigma, [1], [2,2]))/log(2)) == chi
cvx_end


