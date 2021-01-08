function cvx_optval = tight_upper_bound(sigma, chi)

n = 4;

cvx_begin
    variable rho(4,4) hermitian semidefinite
    maximize quantum_cond_entr(rho, [2 2]) / log(2)
    
    subject to 
        trace(rho) == 1
        % (- trace_logm(sigma, rho) + trace_logm(TrX(sigma, [1], [2,2]), TrX(rho, [1], [2,2]))) / log(2) == chi
        % - trace_logm(sigma, rho) / log(2) == chi
        % trace_logm(rho, sigma) / log(2) >= 0.2
        % - trace(rho * logm(sigma)/log(2)) >=0.1 
        
        -trace(rho * logm(sigma)/log(2)) + trace(TrX(rho, [1], [2,2]) * logm(TrX(sigma, [1], [2,2]))/log(2)) == chi

cvx_end


