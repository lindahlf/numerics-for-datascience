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
for tau=1:length(node_index_vector)
    Itau=node_index_vector{tau};
    if (length(Itau)== 1 )
        b(Itau)=b(Itau)+A(Itau,?? )*x(??);
    else
        child1=child_index(tau).child1_number;
        child2=child_index(tau).child2_number;
        Isigma1=node_index_vector{child1};
        Isigma2=node_index_vector{child2};


        block1_add=A(Isigma1,??)*x(Isigma2));
        block2_add=??*x(??);
        b(Itau)=b(Itau)+[block1_add;block2_add];
    end
end
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
