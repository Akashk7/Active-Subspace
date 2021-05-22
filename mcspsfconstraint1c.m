function [Psf,Pf,Re,b]= mcspsfconstraint1c(xmin,lbx,ubx,sd,Nmcs,Pftarget)

[M,~]=size(xmin);

ind = Pftarget*Nmcs;
t0 = xmin;
sd1 = sd;
lbx1 = lbx;
ubx1 = ubx;

%% Calculating Probability

Psf = zeros(M,1);
Pf = zeros(M,1);
Re = zeros(M,1);
b = zeros(M,1);
for k=1:M
x = zeros(Nmcs,2);
mu1 = zeros(1,2);
for i=1:2
mu1(1,i)=t0(k,i);
pd = makedist('normal',mu1(1,i),sd1(1,i));
T = truncate(pd,lbx1(1,i),ubx1(1,i));
x(:,i)=random(T,Nmcs,1);
end

R = constraint1c(x);
Pfbs = R;

[Pfbs]=sortrows(Pfbs,'ascend');

Psf(k,1) = (Pfbs(ind,1));

Pf(k,1) = mean(R<1);
Re(k,1) = 1-Pf(k,1);
b(k,1) = -norminv(Pf(k,1));

end

end