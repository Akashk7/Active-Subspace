function [Psf,Pf,Re,b]= mcspsfconstraint2c(xmin,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget)

[M,~]=size(xmin);

ind = Pftarget*Nmcs;

%% Calculating Probability

mu1=zeros(1,3);
sd1 = [sd(1) sd(2) sd(3)];
lbx1 =[lbx(1) lbx(2) lbx(3)];
ubx1 = [ubx(1) ubx(2) ubx(3)];
Psf = zeros(M,1);
t0 = [xmin(:,1) xmin(:,2) xmin(:,3)];
Pf = zeros(M,1);
Re = zeros(M,1);
b = zeros(M,1);
for k=1:M
x = zeros(Nmcs,5);
for i=1:3
mu1(1,i)=t0(k,i);
pd = makedist('normal',mu1(1,i),sd1(1,i));
T = truncate(pd,lbx1(1,i),ubx1(1,i));
x(:,i)=random(T,Nmcs,1);
end

lbp = mup-3*sdp;
ubp = mup+3*sdp;

pdp1 = makedist('normal',mup(1),sdp(1));
Tp1 = truncate(pdp1,lbp(1),ubp(1));

pdp3 = makedist('normal',mup(3),sdp(3));
Tp3 = truncate(pdp3,lbp(3),ubp(3));

x(:,4)= random(Tp1,Nmcs,1);
x(:,5)= random(Tp3,Nmcs,1);

R = constraint2c(x);

Pfbs= R;

[Pfbs]=sortrows(Pfbs,'ascend');

Psf(k,1) = mean(Pfbs(ind,1));
Pf(k,1) = mean(R<1);
Re(k,1) = 1-Pf(k,1);
b(k,1) = -norminv(Pf(k,1));
end
end