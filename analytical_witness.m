function funs = analytical_witness
    funs.cvenn_werner_like = @cvenn_werner_like
    funs.state_on_boundary = @state_on_boundary;
    funs.grad = @grad;
    funs.witness = @witness;
    funs.werner_like_state = @werner_like_state;
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

function g = grad(sigma, d)
    % need to take care of log 0.
    I = eye(d)
    sigma_b = TrX(sigma, 1, [d,d]);
    g = transpose(logm(sigma) - kron(I, logm(sigma_b)));
end
    
function w = witness(rho, d)
    sigma = state_on_boundary(rho, d)
    g = grad(sigma, d);
    I = eye(d*d);
    w = (trace(sigma*g)*I - g)/(sqrt(trace(g)^2));
end

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
