function [c,ceq] = nonlconinactive_c01_1(x,t1,W1,W21,W22,t21,t22)
c=[];
ceq=[x*W1-t1;x*W21-t21;x*W22-t22];
end
