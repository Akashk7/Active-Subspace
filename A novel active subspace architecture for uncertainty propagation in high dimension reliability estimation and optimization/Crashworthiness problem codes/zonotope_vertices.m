function[Y,X]= zonotope_vertices(W1,Nsamples,maxcount)
%  %"""Compute the vertices of the zonotope.
%     Parameters
%     ----------
%     W1 : ndarray 
%         m-by-n matrix that contains the eigenvector bases of the n-dimensional 
%         active subspace
%     Nsamples : int, optional
%         number of samples per iteration to check (default 1e4)
%     maxcount : int, optional
%         maximum number of iterations (default 1e5)
%     Returns
%     -------
%     Y : ndarray 
%         nzv-by-n matrix that contains the zonotope vertices
%     X : ndarray 
%         nzv-by-m matrix that contains the corners of the m-dimensional hypercube
%         that map to the zonotope vertices
%     %""
[m,n] = size(W1);
totalverts = nzv(m,n);

Z = normrnd(0,1,Nsamples,n);
X  = unique(sign(Z*W1'),'rows');
X = unique(vertcat(X,-X),'rows');
N = size(X,1);

count = 0;
while N<totalverts
    Z = normrnd(0,1,Nsamples,n);
    X0 = unique(sign(Z*W1'),'rows');
    X0 = unique(vertcat(X0,-X0),'rows');
    X = unique(vertcat(X,X0),'rows');
    N = size(X,1);
    count = count+1;
    if count>maxcount
        break
    end
end
numverts = size(X,1);
Y = X*W1;
Y = reshape(Y,numverts,n);
X = reshape(X,numverts,m);
end