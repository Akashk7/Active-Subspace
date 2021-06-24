function [c,ceq]=nonlcon_approxmoments(x,srgt,srgt_evaluate,beta)
c=[];
ceq = srgt_evaluate(x, srgt)-beta;
end