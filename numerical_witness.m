function ans = numerical_witness(rho, d)

    if d==2
        s = closest_cvenn_state(rho)
    else 
        s = closest_cvenn_state_3D(rho)
    end
    I = eye(d^2)
    %ans = trace((I - rho/trace(rho*s)) * x)
    ans = (trace(s*rho - s*s)*I + s - rho)/norm(rho - s, 'fro')