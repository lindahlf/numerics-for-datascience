function b = matvec(x,A,child_index,node_index_vector)
b=zeros(n,1);

for tau=1:length(node_index_vector)
    Itau=node_index_vector{tau};
    if (length(Itau)== 1 )
        b(Itau)=b(Itau)+A(Itau,Itau )*x(Itau);
    else
        child1=child_index(tau).child1_number;
        child2=child_index(tau).child2_number;
        Isigma1=node_index_vector{child1};
        Isigma2=node_index_vector{child2};

        block1_add=A(Isigma1,Isigma2)*x(Isigma2);
        block2_add=A(Isigma2,Isigma1)*x(Isigma1);
        b(Itau)=b(Itau)+[block1_add;block2_add];
    end
end
end

