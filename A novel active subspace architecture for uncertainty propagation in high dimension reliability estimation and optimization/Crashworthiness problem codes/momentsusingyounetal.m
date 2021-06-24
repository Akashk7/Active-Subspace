function [y_mean,y_std] = momentsusingyounetal(y)
M = size(y,1);
P  = zeros(M,1);
beta = zeros(M,1);
y1 = sortrows(y,'ascend');
for i=1:M
    P(i,1) = (i-0)/(M+1);
    beta(i,1) = norminv(P(i,1));
end

figure;
plot(y1,beta,'ro--');
xlabel('PSF');
ylabel('beta');

beta_1 = -sqrt(3);
beta_2 = 0;
beta_3 = sqrt(3);

a0 = [1 1 1 1];
lba = [0 0 0 0];
uba = [2 2 2 2];
func = @(a)(norm([ones(M,1) beta beta.^2 beta.^3]*[a(1);a(2);a(3);a(4)]-y1));
[amin,fval,exitflag]=fmincon(func,a0,[],[],[],[],lba,uba,@(a)nonlcona(a));

y1_1 = [1 beta_1 beta_1.^2 beta_1.^3]*[amin(1);amin(2);amin(3);amin(4)];
y1_2 = [1 beta_2 beta_2.^2 beta_2.^3]*[amin(1);amin(2);amin(3);amin(4)];
y1_3 = [1 beta_3 beta_3.^2 beta_3.^3]*[amin(1);amin(2);amin(3);amin(4)];


y_mean = (1/6*(y1_1+y1_3))+(4/6*(y1_2));
y_std = sqrt((1/6*((y1_1-y_mean).^2))+(1/6*((y1_3-y_mean).^2)));

