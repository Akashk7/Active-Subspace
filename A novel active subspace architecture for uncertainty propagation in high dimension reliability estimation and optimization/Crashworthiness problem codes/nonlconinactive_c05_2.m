function [c,ceq] = nonlconinactive_c05_2(x,t1,W1,W21,t21)
c=[];
ceq=[x*W1-t1;x*W21-t21];
end
