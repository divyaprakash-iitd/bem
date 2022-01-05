function [x,y,phi] = calculate_domain(bem,N)
% Calculates the unknowns at N points distributed uniformly inside the
% domain
    r = linspace(0,0.9,N);
    theta = linspace(0,2*pi,N);
    [R, THETA] = meshgrid(r,theta);
    x = R.*cos(THETA);
    y = R.*sin(THETA);
    phi = zeros(N,N);
    % Calculate the solution
    for iN = 1:N
        for jN = 1:N
            xi  = x(iN,jN);
            eta = y(iN,jN);
            phi(iN,jN) = sol_in(bem,xi,eta);
        end
    end
end
