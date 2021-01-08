function cvx_optval = tight_upper_bound(sigma, chi)

n = 4;

cvx_begin
    variable rho(4,4) hermitian semidefinite
    maximize quantum_cond_entr(rho, [2 2]) / log(2)
    t = class(sigma)
    
    subject to 
        trace(rho) == 1
        % (- trace_logm(sigma, rho) + trace_logm(TrX(sigma, [1], [2,2]), TrX(rho, [1], [2,2]))) / log(2) == chi
        - trace_logm(cvx_constant(sigma), rho) / log(2) == chi
cvx_end


