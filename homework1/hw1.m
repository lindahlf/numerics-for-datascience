%% Part 2b
clear all; close all; clc

V = [1 2 2003 2005; 2 2 2002 2004; 3 2 2001 2003; 4 7 7005 7012];

[m n] = size(V);
Q = []; R = []; A = V; p = min(m,n); A_error = [];

for i = 1:p
    q = A(:,i)/norm(A(:,i));
    rt = q'*A;
    Q = [Q q];
    R = [R; rt];
    A = A - q*rt;
    A_error = [A_error norm(A)];
end

plot([1:1:p],A_error)

%% Part 2c
clear all; close all; clc

V = load_mat_hw1(1000,100);

[m n] = size(V);
Q = []; R = []; A = V; p = min(m,n); A_error = [];

for i = 1:p
    q = A(:,i)/norm(A(:,i));
    rt = q'*A;
    Q = [Q q];
    R = [R; rt];
    A = A - q*rt;
    A_error = [A_error norm(A)];
end

plot([1:1:p],A_error)

%% Part 2d
clear all; close all; clc

V = load_mat_hw1(1000,100);

[m n] = size(V);
Q = []; R = []; A = V; p = min(m,n); A_error = [];

for j = 1:p
    [~, i] = max(vecnorm(A));
    q = A(:,i)/norm(A(:,i));
    rt = q'*A;
    Q = [Q q];
    R = [R; rt];
    A = A - q*rt;
    A_error = [A_error norm(A)];
end

plot([1:1:p],A_error)

%% Part 4
clear all; close all; clc

W = load_mat_hw1(1000,100);

[m n] = size(W);
Q = []; R = []; A = W; p = min(m,n); 

k = 1;

while norm(A) > 10e-10
    [~, i] = max(vecnorm(A));
    q = A(:,i)/norm(A(:,i));
    rt = q'*A;
    Q = [Q q];
    R = [R; rt];
    A = A - q*rt;
    k = k + 1;
end

[U_hat, S, V] = svd(R);
U = Q*U_hat;


