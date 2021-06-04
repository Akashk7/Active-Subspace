function N1 = nzv_2(m,n)    
% """Number of zonotope vertices.
%     
%     Compute the number of zonotope vertices for a linear map from R^m to R^n.
%     
%     Parameters
%     ----------
%     m : int
%         the dimension of the hypercube
%     n : int
%         the dimension of the low-dimesional subspace
%         
%     Returns
%     -------
%     N : int 
%         the number of vertices defining the zonotope
%     """
N = 0;
for i=1:n-1
    N = N + nchoosek(m-1,i);
end
N1 = 2*N;
