function [zmin,zmax] = inactivebnd(t,W1,W2)
[m,~]=size(W1);
flb = 1;
fub = -1;
Al = [-W2;W2];
y = t';
bl = [((ones(m,1))+W1*y);((ones(m,1))-W1*y)];
func = @(z)(flb*z);
splb= 1;
zmin=fmincon(func,splb,Al,bl);

func_max = @(z)(fub*z);
spub=1;
zmax=fmincon(func_max,spub,Al,bl);

end
