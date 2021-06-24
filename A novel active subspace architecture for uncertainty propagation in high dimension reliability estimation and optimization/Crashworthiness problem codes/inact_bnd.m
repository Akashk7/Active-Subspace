function [zmin,zmax] = inact_bnd(t1,W1,W2)
[M,m] = size(W2);

[zmin,~,exitflag_zmin] = linprog(ones(1,m),[-W2;W2],[W1*t1'+ones(M,1);ones(M,1)-W1*t1']);
[zmax,~,exitflag_zmax] = linprog(-1*ones(1,m),[-W2;W2],[W1*t1'+ones(M,1);ones(M,1)-W1*t1']);
end


