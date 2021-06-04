function [Psf,Pf,Re,b]= mcspsfconstraint4c(xmin,lbx,ubx,sd,Nmcs,Pftarget)

[M,~]=size(xmin);

ind = Pftarget*Nmcs;

%% Calculating Probability
t0 = xmin(:,1:17);
sd1 = sd(:,1:17);
lbx1 = lbx(:,1:17);
ubx1 = ubx(:,1:17);
Psf = zeros(M,1);
Pf = zeros(M,1);
Re = zeros(M,1);
b = zeros(M,1);

for k=1:M
x = zeros(Nmcs,17);
mu1 = zeros(1,17);

for i=1:17
    if i<16
        mu1(1,i)=t0(k,i);
        pd = makedist('normal',mu1(1,i),sd1(1,i));
        T = truncate(pd,lbx1(1,i),ubx1(1,i));
        x(:,i)=random(T,Nmcs,1);
    elseif i>15
        mu1(1,i)=xmin(k,i);
        pd1 = makedist('Lognormal',mu1(1,i),sd1(1,i));
        T1 = truncate(pd1,lbx1(1,i),ubx1(1,i));
        x(:,i)=random(T1,Nmcs,1);
    end
end

R = constraint4c(x);

Pfbs= R;
[Pfbs]=sortrows(Pfbs,'ascend');

Psf(k,1) = (Pfbs(ind,1));
Pf(k,1) = mean(R<1);
Re(k,1) = 1-Pf(k,1);
b(k,1) = -norminv(Pf(k,1));
end

end