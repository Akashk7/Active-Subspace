function [t2_min, t2_max]= inactbnd(t1,W1,W2)
t2_min  = (-1-t1*W1')*pinv(W2');
t2_max = (1-t1*W1')*pinv(W2');
end