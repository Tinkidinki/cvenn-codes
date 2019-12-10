function funs = analytical_witness
    funs.cvenn_werner_like = @cvenn_werner_like
    funs.state_on_boundary = @state_on_boundary;
    % funs.grad = @grad;
    funs.witness = @witness;
    funs.werner_like_state = @werner_like_state;
    funs.test = @testfun;
    funs.non_optimal_witness = @non_optimal_witness;
    funs.witness_performance_on_werner = @witness_performance_on_werner;
    funs.witness_performance_on_isotropic = @witness_performance_on_isotropic;
end

% function g = grad(sigma)
% % need to take care of log 0.
%     I = [1 0; 0 1];
%     sigma_b = TrX(sigma, 1, [2,2]);
%     g = transpose(logm(sigma) - kron(I, logm(sigma_b)));
% end

% function w = witness(sigma)
%     g = grad(sigma);
%     I = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
%     w = (trace(sigma*g)*I - g)/(sqrt(trace(g)^2));
% end

function c = cvenn_werner_like(p, rho, d)
    I = eye(d*d)
    sigma = p*rho + (1-p)*I/d^2
    sigma_b = TrX(sigma, 1, [d,d])
    c = trace(sigma*logm(sigma)) - trace(sigma_b*logm(sigma_b))
end

function w = witness(rho, d)
    % need to take care of log 0.
    sigma = state_on_boundary(rho, d)
    I = eye(d)
    sigma_b = TrX(sigma, 1, [d,d]);
    w = -logm(sigma)/log(2) + kron(I, logm(sigma_b)/log(2));
end
    
% function w = witness(rho, d)
%     sigma = state_on_boundary(rho, d)
%     g = grad(sigma, d);
%     I = eye(d*d);
%     w = (trace(sigma*g)*I - g)/(sqrt(trace(g)^2));
% end

function sigma = state_on_boundary(rho, d)
    p = fsolve(@(p)cvenn_werner_like(p, rho, d), 0.1)
    % sigma = p
    sigma = p*rho + (1-p)*eye(d^2)/d^2
end

function rho = werner_like_state(p,d)
    phi_plus_vector = zeros(d^2,1)
    I_d = eye(d)
    for i = 1:d 
        phi_plus_vector = phi_plus_vector + kron(I_d(:,i),I_d(:,i))
    end
    phi_plus_matrix = kron(phi_plus_vector, phi_plus_vector')
    phi_plus_matrix = phi_plus_matrix/d
    rho = p*phi_plus_matrix + (1-p)*eye(d^2)/d^2
end

% No, need to change, has been written correctly.
function w = non_optimal_witness(rho_s, d)
    I = eye(d)
    rho_sb = TrX(rho_s, 1, [d,d]);
    w = -logm(rho_s)/log(2) + kron(I, logm(rho_sb)/log(2));
end


% function ans = testfun(d)
%     ans = 0
%     for i = 1:100
%         rho = randRho(d*d)
%         if quantum_cond_entr(rho, [d d]) < 0:
%             disp(i);
%             w = witness(rho, d)
%             if check_witness(w, rho) >= 0:
%                 disp("uh oh")
%                 ans = -1

function wp = witness_performance_on_werner(w, d)
    prec = 100; % set precision, the higher this number is, the better - will take longer to compute though.
    p = zeros(1, prec);
    cond_entr = zeros(1, prec);  %The quantum conditional entropy
    wit_res = zeros(1, prec);  %The result of the witness
    for i = 1:prec
        p(i) = i/prec;
        rho = werner_like_state(p(i), d);
        cond_entr(i) = quantum_cond_entr2(rho, d);
        wit_res(i) = check_witness(w, rho);
    end
    
    a1 = plot(p, cond_entr); M1 = "Conditional Entropy"
    hold on;
    a2 = plot(p, wit_res); M2 = "Value of $Tr(W \rho)$"
    leg = legend([a1,a2], [M1, M2]);
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
    
    a1 = plot(p, cond_entr); M1 = "Conditional Entropy"
    hold on;
    a2 = plot(p, wit_res); M2 = "Value of $Tr(W \rho)$"
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

