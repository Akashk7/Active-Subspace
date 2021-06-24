function n = active_subspace_dimension(lambda)
M=size(lambda,1);
difference = zeros(M-1,1);
for i=1:M-1
    difference(i,1) = lambda(i,1)-lambda(i+1,1);
end
[~,I] = sort(difference,'descend');
n = I(1);
end

