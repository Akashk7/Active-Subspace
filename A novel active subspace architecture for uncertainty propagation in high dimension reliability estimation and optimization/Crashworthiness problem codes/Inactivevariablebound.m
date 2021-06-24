function [lbox,ubox] = Inactivevariablebound(t,W1,W2)
[M,~]=size(t);
[m,n]=size(W1);
for j= 1:M
    y = t(j,1:n)';
    Al = [-W2;W2];
    bl = [1+(W1*y);1-(W1*y)];
    for i=1:m-n    
    clb= zeros(m-n,1);
    clb(i,1) = 1;
    [zlb] = linprog(clb,Al,bl);
    lbox(1,i) =zlb(i,1);
    cub = zeros(m-n,1);
    cub(i,1)= -1;
    zup = linprog(cub,Al,bl);
    ubox (1,i)= zup(i,1);
    end
end


end
