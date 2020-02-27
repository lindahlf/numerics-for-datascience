function y = fftx(x)

n = length(x);
omega = exp(-2*1i*pi/n);

if rem(n,2) == 0
    k = (0:n/2-1)'; w = omega.^k;
    u = fftx(x(1:2:n-1));
    v = w.*fftx(x(2:2:n));
    y = [u + v; u - v];
else
    j = 0:n-1; k = j';
    F = omega.^(k*j);
    y = F*x;
end

end

