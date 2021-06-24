function [c,ceq] = nonlconinactive_c04_3(x,t1,W1,W21,W22,W23,W24,t21,t22,t23,t24)
c=[];
ceq=[x*W1-t1;x*W21-t21;x*W22-t22;x*W23-t23;x*W24-t24];
end
