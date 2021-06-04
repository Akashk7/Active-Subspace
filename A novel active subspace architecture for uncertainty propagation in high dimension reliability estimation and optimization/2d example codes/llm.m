function [b,C,W,Ev,Average_R2_pred]=llm(xN,qoi,Minda,p)
M=Minda;
out1 = randperm(M);
ind = out1(1:M);
ind = ind';
xM=xN(ind,:);
N=size(xN,1);
dist = zeros(N,M);
    for i=1:M
        diff(1:N,:) = repmat(xM(i,:),N,1) - xN;
        d1 = diff.^2;
        d2 = d1';
        d3 = sum(d2);
        dist(1:N,i) = d3';
    end
[D,I]=sort(dist,'ascend');
ind = (D~=0);
I1 = I(ind);
I1 = reshape(I1,N-1,M);
for i=1:M
    X = sqrt(1./N)*[ones(p,1) xN(I1(1:p,i),:)];
    Y = sqrt(1./N)*qoi(I1(1:p,i),1); 
    B = pinv(X)*Y;
    b(i,:) = B(2:end,1)';
    [srgt_PRS, PRESSRMS_PRS, eXV_PRS, srgtOPT_PRS, Y_hat_PRS, pred_var_PRS, R2_pred_PRS(i,:)] = metamodel_PRS(X,Y);
end

Average_R2_pred = mean(R2_pred_PRS);

%% Calculating C Matrix
C = 1/M*(b'*b);

%% Eigendecomposition of C

[EigenvectorsC,~]=eig(C);
W = EigenvectorsC(:,end:-1:1);
[~,Ev,~]=svd(sqrt(1./M).*b');
Ev=(Ev).^2;


