function [c,ceq] = nonlconinactive_c04_2(x,t1,W1,W22,W23,t22,t23)
c=[];
ceq=[x*W1(:,1)-t1(1);x*W1(:,2)-t1(2);x*W1(:,3)-t1(3);x*W22-t22;x*W23-t23];
end
