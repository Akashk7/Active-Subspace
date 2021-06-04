function [c,ceq] = nonlconinactive(x,t1,W1,W2,t2)
[~,n1]=size(W1);
[~,n2]=size(W2);
m=n1+n2;
W = [W1 W2];
t = [t1 t2];
c=[];
ceq=zeros(m,1);
for i=1:m
ceq(i,1) = x*W(:,i)-t(i);
end
end
