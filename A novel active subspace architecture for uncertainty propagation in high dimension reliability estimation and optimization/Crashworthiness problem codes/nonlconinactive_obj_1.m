function [c,ceq] = nonlconinactive_obj_1(x,t1,W1,W22,W23,W24,W25,t22,t23,t24,t25)
c=[];
ceq=[x*W1-t1;x*W22-t22;x*W23-t23;x*W24-t24;x*W25-t25];
end
