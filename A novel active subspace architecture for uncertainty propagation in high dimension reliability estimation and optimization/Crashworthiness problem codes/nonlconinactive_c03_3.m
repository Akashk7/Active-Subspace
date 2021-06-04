function [c,ceq] = nonlconinactive_c03_3(x,t1,W1,W3,t3)
c=[];
ceq=[x*W1(:,1)-t1(1);x*W1(:,2)-t1(2);x*W1(:,3)-t1(3);x*W3-t3];
end
