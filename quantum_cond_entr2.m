function cond_entr = quantum_cond_entr2(rho, d)
    rho_b = TrX(rho, 1, [d d]);

    Term1 = 0;
    Term2 = 0;

    for i = eig(rho)'
        if (i~=0)
            Term1 = Term1 + (i*log2(i))
        end
    end

    for i = eig(rho_b)'
        if(i~=0)
            Term2 = Term2 + (i*log2(i))
        end
    end

    cond_entr = - Term1 + Term2
end
