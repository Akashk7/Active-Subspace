function [c,ceq] = nonlconinactive_c03_4(x,t1,W1,W21,W22,W23,t21,t22,t23)
c=[];
ceq=[x*W1-t1;x*W21-t21;x*W22-t22;x*W23-t23];
end
