function [c,ceq] = nonlconinactive_c08_1(x,t1,W1,W21,t2)
c=[];
ceq=[x*W1(:,1)-t1(1);x*W1(:,2)-t1(2);x*W21-t2];
end
