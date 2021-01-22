function cvx_optval = tight_upper_bound_witness(witness, chi)

n = 4;

cvx_begin
    variable rho(4,4) hermitian semidefinite
    maximize quantum_cond_entr(rho, [2 2]) / log(2)
    
    subject to 
        trace(rho) == 1
        trace(witness * rho) == chi
cvx_end


