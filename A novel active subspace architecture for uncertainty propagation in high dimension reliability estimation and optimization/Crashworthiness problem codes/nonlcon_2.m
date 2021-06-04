function [c,ceq]=nonlcon_2(x,t1,t2_comb,W1,W2)
c=[];
W = [W1 W2];
t = [t1 t2_comb];
ceq = x*W-t;
end