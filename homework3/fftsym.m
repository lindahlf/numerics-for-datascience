function y = fftsym(x, symmetry)


n = length(x);
omega = exp(-2*1i*pi/n);

if rem(n,2) == 0
    k = (0:n/2-1)'; w = omega.^k;
    u = fftsym(x(1:2:n-1),0);
    if symmetry == 1        
        v = w.*u;
    else
        v = w.*fftsym(x(2:2:n),0);
    end
    y = [u + v; u - v];
else
    j = 0:n-1; k = j';
    F = omega.^(k*j);
    y = F*x;
end

end

% Difference for this function should be that we use computations
% 