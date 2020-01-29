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

k = 0;

while norm(A) > 1e-10
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

%% Part 6d
clear all; close all; clc

A = zeros(690*928*3,27);

for k = 1:27
    filename = sprintf("testbild_snapshots_%04d.png",k);
    frame_k = imread(filename);
    imfloat = double(frame_k);
    n1 = size(imfloat,1); n2 = size(imfloat,2); n3 = 3;
    n=n1*n2*n3;
    
    v = reshape(imfloat,n,1);
    
    A(:,k) = v;
end

tic
[U, S, V] = svd(A,0);
toc

u = A(:,1);
v = ones(27,1);

test = A - u*v';

%% Part 6e 
% Compute an approximate SVD of the matrix A
clear all; close all; clc

% Generate A 
A = zeros(690*928*3,27);

for k = 1:27
    filename = sprintf("testbild_snapshots_%04d.png",k);
    frame_k = imread(filename);
    imfloat = double(frame_k);
    n1 = size(imfloat,1); n2 = size(imfloat,2); n3 = 3;
    n=n1*n2*n3;
    
    v = reshape(imfloat,n,1);
    
    A(:,k) = v;
end

tic
% Start of SVD approx. 
W = A; 

[m n] = size(W);
Q = []; R = []; A = W; p = min(m,n); 

k = 0;

for k = 0:2
    [M, i] = max(vecnorm(A));
    q = A(:,i)/M;
    rt = q'*A;
    Q = [Q q];
    R = [R; rt];
    A = A - q*rt;
    %k = k + 1;
end

[U_hat, S, V] = svd(R);
U = Q*U_hat;
toc

%% Part 7
clear all; close all; clc

k = 1;
filename = sprintf("roundabout_snapshots_%04d.png",k);
test_frame = imread(filename);
imfloat = double(test_frame);
n1 = size(imfloat,1); n2 = size(imfloat,2); n3 = 3;
n = n1*n2*n3;

A = zeros(n,56);

for k = 1:56
    filename = sprintf("roundabout_snapshots_%04d.png",k);
    frame_k = imread(filename);
    imfloat = double(frame_k);
    
    n1 = size(imfloat,1); n2 = size(imfloat,2); n3 = 3;
    n = n1*n2*n3;
    
    v = reshape(imfloat,n,1);
    
    A(:,k) = v;
end

tic
% Start of SVD approx. 
W = A; 

[m n] = size(W);
Q = []; R = []; A = W; p = min(m,n); 

A_error = [];

for k = 0:10
    [M, i] = max(vecnorm(A));
    q = A(:,i)/M;
    rt = q'*A;
    Q = [Q q];
    R = [R; rt];
    A = A - q*rt;
    
end

[U_hat, S, V] = svd(R);
U = Q*U_hat;
toc
%% computing b1

u = U(:,1);
v = V(:,1);
s = S(1,1);

b1 = s*u*v(1);


imb = reshape(b1, n1, n2, n3);
imshow(uint8(imb))


%% Part 8
clear all; close all; clc


filename = sprintf("roundabout_snapshots_%04d.png",1);
test_frame = imread(filename);
imfloat = double(test_frame);
n1 = size(imfloat,1); n2 = size(imfloat,2); n3 = 3;
n = n1*n2*n3;

A = zeros(n,56); %56 is the number of frames

for k = 1:56
    filename = sprintf("roundabout_snapshots_%04d.png",k);
    frame_k = imread(filename);
    imfloat = double(frame_k);
    
    n1 = size(imfloat,1); n2 = size(imfloat,2); n3 = 3;
    n = n1*n2*n3;
    
    v = reshape(imfloat,n,1);
    
    A(:,k) = v;
end

tic
k = 1; p = 4; 
G = random('normal',0,1,[56,k+p]);

Y = A*G; Q = orth(Y);
B = Q'*A;
[U_hat, D, V] = svd(B);
U = Q*U_hat;
toc

% Computing b1 (again)
u = U(:,1);
v = V(:,1);
d = D(1,1);

b1 = d*u*v(1);
imb = reshape(b1, n1, n2, n3);
imshow(uint8(imb))



