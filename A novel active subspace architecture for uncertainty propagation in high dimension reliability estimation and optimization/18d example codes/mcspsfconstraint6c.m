function [Psf,Pf,Re,b]= mcspsfconstraint6c(xmin,lbx,ubx,sd,Nmcs,Pftarget)

[M,~]=size(xmin);

ind = Pftarget*Nmcs;

%% Calculating Probability
t0 = xmin;
sd1 = sd;
lbx1 = lbx;
ubx1 = ubx;
Psf = zeros(M,1);
Pf = zeros(M,1);
Re = zeros(M,1);
b = zeros(M,1);
for k=1:M
x = zeros(Nmcs,18);
mu1 = zeros(1,18);

for i=1:15
mu1(1,i)=t0(k,i);
pd = makedist('normal',mu1(1,i),sd1(1,i));
T = truncate(pd,lbx1(1,i),ubx1(1,i));
x(:,i)=random(T,Nmcs,1);
end

for i=16:18
mu1(1,i)=xmin(k,i);
pd1 = makedist('Lognormal',mu1(1,i),sd1(1,i));
T1 = truncate(pd1,lbx1(1,i),ubx1(1,i));
x(:,i)=random(T1,Nmcs,1);
end

R = constraint6c(x);

Pfbs = R;

[Pfbs]=sortrows(Pfbs,'ascend');

Psf(k,1) = mean(Pfbs(ind,1));
Pf(k,1) = mean(R<1);
Re(k,1) = 1-Pf(k,1);
b(k,1) = -norminv(Pf(k,1));
end

end