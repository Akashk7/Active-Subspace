function [c,ceq] = nonlconinactive_c02_1(x,t1,W1,W22,t2)
c=[];
ceq=[x*W1-t1;x*W22-t2];
end
