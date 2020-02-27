function [tree,child_vec,node_index_vec,number]=build_bisection_tree(I,child_vec,node_index_vec)
% This builds a bisection tree of the indices in I.
% The function returns:
%    child_vec: a vector of structs, such that child_vec(k).child1
%               and child_vec(k).child2 returns the numbers for the
%               children of node k
%    node_index_vec: a cell array containing indices. The indices
%                    for node k are given by node_index_vec(k).
%

    s0=min(I);
    s1=max(I);
    smean=floor(median(I));
    if (s1>s0)
        [child1,child_vec,node_index_vec,number1]=build_bisection_tree(s0:smean,child_vec,node_index_vec);
        [child2,child_vec,node_index_vec,number2]=build_bisection_tree(smean+1:s1,child_vec,node_index_vec);
    else
        child1=NaN;
        child2=NaN;
        number1=NaN;
        number2=NaN;
    end
    node_index_vec={node_index_vec{:},I};


    child_vec=[child_vec(:); struct('child1_number',number1,'child2_number',number2)];
    tree=struct('index_vector',I,'child1',child1,'child2',child2);
    number=length(node_index_vec);
