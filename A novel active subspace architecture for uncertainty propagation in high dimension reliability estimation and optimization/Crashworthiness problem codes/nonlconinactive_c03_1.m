function [c,ceq] = nonlconinactive_c03_1(x,t1,W1,W3)
c=[];
ceq=[x*W1-t1;x*W3];
end
