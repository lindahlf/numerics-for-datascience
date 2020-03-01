% Compute a semiseparable matrix of size n times n
n=500;
randn('seed',0);
v=randn(n,1);
D=diag(randn(n,1));
A=v*v'+D;

% Build a bisection tree
[~,child_index,node_index_vector]=build_bisection_tree(1:n,[],{});
disp('Computed tree');


b=zeros(n,1);
x=randn(n,1);
tic; % Start timing here since the building of the tree is not optimized

tt1=toc;

disp('completed structured matvec');

% Compare with naive matrix vector multiplication
tic
b2=zeros(n,1);
for i=1:n
    for j=1:n
        b2(i)=b2(i)+A(i,j)*x(j);
    end
end
tt2=toc;

disp('completed naive matvec');

norm(b-b2) % should be zero

tt1
tt2
