function [c,ceq] = nonlconinactive_c03_2(x,t1,W1,W21,W23)
c=[];
ceq=[x*W1-t1;[x*W21 x*W23]];
end
