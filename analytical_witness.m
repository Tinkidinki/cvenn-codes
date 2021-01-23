function funs = analytical_witness
    funs.cve_werner_like = @cve_werner_like
    funs.state_on_boundary = @state_on_boundary;
    funs.witness = @witness;
    funs.werner_like_state = @werner_like_state;
    funs.test = @testfun;
    funs.non_optimal_witness = @non_optimal_witness;
    funs.witness_performance_on_werner = @witness_performance_on_werner;
    funs.witness_performance_on_isotropic = @witness_performance_on_isotropic;
    funs.cve_vs_upperbound_werner = @cve_vs_upperbound_werner;
end

% All logarithms have been taken to the base 2 for calculation of entropy. 
% That is, we take logm/log2


function c = cve_werner_like(p, rho, d)
% Generates a state of the form p*rho + (1-p)*I/d^2
% The conditional von Neumann entropy of the generated state is returned. 
    I = eye(d*d)
    sigma = p*rho + (1-p)*I/d^2
    c = quantum_cond_entr2(sigma, d)
end

function w = witness(rho, d)
    % Finds the "state on boundary" for the given state
    % Calculates the witness using that state as the parameter
    sigma = state_on_boundary(rho, d)
    I = eye(d)
    sigma_b = TrX(sigma, 1, [d,d]);
    w = -logm(sigma)/log(2) + kron(I, logm(sigma_b)/log(2));
end
    

function sigma = state_on_boundary(rho, d)
% Returns a state that lies on the line p*rho + (1-p)*I/d^2 and on the boundary of CVENN
    p = fsolve(@(p)cve_werner_like(p, rho, d), 0.1)
    % sigma = p
    sigma = p*rho + (1-p)*eye(d^2)/d^2
end

function rho = werner_like_state(p,d)
% Returns a state of the form (p*|phi+><phi+| + (1-p)I/d^2) where d is the number of qubits the state is defined over, and phi+ is the Bell state generalised to higher dimensions (defined in the paper)
    phi_plus_vector = zeros(d^2,1)
    I_d = eye(d)
    for i = 1:d 
        phi_plus_vector = phi_plus_vector + kron(I_d(:,i),I_d(:,i))
    end
    phi_plus_matrix = kron(transpose(phi_plus_vector), phi_plus_vector)
    phi_plus_matrix = phi_plus_matrix/d
    rho = p*phi_plus_matrix + (1-p)*eye(d^2)/d^2
end

function w = non_optimal_witness(rho_s, d)
% Prepares a witness by substituting the state rho_s in the witness formula. This witness may not be optimal. 
    I = eye(d)
    rho_sb = TrX(rho_s, 1, [d,d]);
    w = -logm(rho_s)/logm(2)+ kron(I, logm(rho_sb)/logm(2));
end

function ub = cve_vs_upperbound_werner(w, d)
% Generates Werner states on teh line p|phi+><phi+| + (1-p)I/d^2
% Finds the expectation value chi of the witness for each state
% Finds the value of upper bound of T_w^chi for each state
% Finds the conditional entropy for each state
% Returns a difference array of upperbound - conditional entropy. 

    prec = 100;
    p = zeros(1, prec);
    cond_entr = zeros(1, prec);
    chi = zeros(1, prec);
    upper_bound = zeros(1, prec);
    difference = zeros(1, prec);

    for i = 1:prec
        p(i) = i/prec;
        rho = werner_like_state(p(i), d);
        cond_entr(i) = quantum_cond_entr(rho, [d d]) / logm(2);
        chi(i) = check_witness(w, rho);
        upper_bound(i) = tight_upper_bound_witness(w, chi(i));
        difference(i) = upper_bound(i) - cond_entr(i);
    end
    disp("Difference between upperbound and conditional entropy")
    disp(difference)
end


function wp = witness_performance_on_werner(w, d)
% Generates Werner states on the line p|phi+><phi+| + (1-p)I/d^2
% Finds the expectation of each state against the witness and plots this. 
% The actual value of conditional entropy is also plotted. 
    prec = 500; % set precision, the higher this number is, the better - will take longer to compute though.
    p = zeros(1, prec);
    cond_entr = zeros(1, prec);  %The quantum conditional entropy
    wit_res = zeros(1, prec);  %The result of the witness
    for i = 1:prec
        p(i) = i/prec;
        rho = werner_like_state(p(i), d);
        cond_entr(i) = quantum_cond_entr2(rho, d);
        wit_res(i) = check_witness(w, rho);
    end
    
    a1 = plot(p, cond_entr, '-'); M1 = "Conditional entropy $\approx$ Upper bound of $T_{W_{\rho wer}}^\chi$"
    a1.LineWidth = 1.5
    hold on;

    a2 = plot(p, wit_res, '--'); M2 = "Value of $\chi$"
    a2.LineWidth = 1.5
    hold on;


    leg = legend([a1,a2], [M1, M2 ]);
    leg.FontSize = 14;
    set(leg, 'Interpreter', 'latex');
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlh = xlabel('$p$: Werner Mixing Parameter', 'interpreter', 'latex');
    xlh.Position(2) = xlh.Position(2)-0.95;
    xlh.Position(1) = xlh.Position(1)-0.3;
    wp = 0;

end

function wi = witness_performance_on_isotropic(w, d)
% Generates Isotropic states (defined in the paper)
% Finds the expectation of each state against the witness and plots this. 
% The actual value of conditional entropy is also plotted. 
    prec = 100; % set precision, the higher this number is, the better - will take longer to compute though.
    p = zeros(1, prec);
    cond_entr = zeros(1, prec);  %The quantum conditional entropy
    wit_res = zeros(1, prec);  %The result of the witness
    for i = 1:prec
        p(i) = -1/(d^2 - 1) + i/prec*(1 + 1/(d^2 - 1));
        rho = werner_like_state(p(i), d);
        cond_entr(i) = quantum_cond_entr2(rho, d);
        wit_res(i) = check_witness(w, rho);
    end
    
    a1 = plot(p, cond_entr, '-'); M1 = "Conditional Entropy"
    a1.LineWidth = 1.5
    hold on;
    a2 = plot(p, wit_res, '--'); M2 = "Value of $Tr(W \rho)$"
    a2.LineWidth = 1.5
    leg = legend([a1,a2], [M1, M2]);
    set(leg, 'Interpreter', 'latex');
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlh = xlabel('$\alpha$: Isotropic Mixing Parameter', 'interpreter', 'latex');
    xlh.Position(2) = xlh.Position(2)-0.95;
    xlh.Position(1) = xlh.Position(1)-0.3;
    wi = 0;

end

