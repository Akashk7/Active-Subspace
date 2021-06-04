function [Psf,Pf,Re,b]= mcspsfconstraint10c(xmin,lbx,ubx,sd,mup,sdp,Nmcs,Pftarget)

[M,~]=size(xmin);

ind = Pftarget*Nmcs;

%% Calculating Probability

t0 = [xmin(:,3) xmin(:,5) xmin(:,6) xmin(:,7)];
mu1 = zeros(1,4);
sd1 = [sd(3) sd(5) sd(6) sd(7)];
lbx1 = [lbx(3) lbx(5) lbx(6) lbx(7)];
ubx1 = [ubx(3) ubx(5) ubx(6) ubx(7)];
Psf = zeros(M,1);
Pf = zeros(M,1);
Re = zeros(M,1);
b = zeros(M,1);

for k=1:M

x = zeros(Nmcs,7);

for i=1:4
mu1(1,i)=t0(k,i);
pd = makedist('normal',mu1(1,i),sd1(1,i));
T = truncate(pd,lbx1(1,i),ubx1(1,i));
x(:,i)=random(T,Nmcs,1);
end

lbp = mup-3*sdp;
ubp = mup+3*sdp;

pdp2 = makedist('normal',mup(2),sdp(2));
Tp2 = truncate(pdp2,lbp(2),ubp(2));

pdp3 = makedist('normal',mup(3),sdp(3));
Tp3 = truncate(pdp3,lbp(3),ubp(3));

pdp4 = makedist('normal',mup(4),sdp(4));
Tp4 = truncate(pdp4,lbp(4),ubp(4));


x(:,5) = random(Tp2,Nmcs,1);
x(:,6) = random(Tp3,Nmcs,1);
x(:,7) = random(Tp4,Nmcs,1);


R = constraint10c(x);
Pfbs = R;

[Pfbs]=sortrows(Pfbs,'ascend');

Psf(k,1) = mean(Pfbs(ind,1));
Pf(k,1) = mean(R<1);
Re(k,1) = 1-Pf(k,1);
b(k,1) = -norminv(Pf(k,1));

end

end