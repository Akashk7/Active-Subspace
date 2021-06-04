function [c,ceq] = nonlconinactive_c06_1(x,t1,W1,W22,t22)
c=[];
ceq=[x*W1(:,1)-t1(1);x*W1(:,2)-t1(2);x*W1(:,3)-t1(3);x*W1(:,4)-t1(4);x*W22-t22];
end
