% script to draw distribution functions for the mean dwell time as a
% function of average mean dwell time and strength

tau=3;
var=100*tau;

w1=1+tau*(tau-1)/var;
w2=(tau-1)*w1;

x=logspace(0,4,1000);
y=x.^-(w1+w2-1).*(x-1).^(w2-1)/beta(w1,w2);

disp(num2str([w1 w2]))

figure(1)
%clf
hold on
plot(log(x),y,'k-')
set(gca,'xscale','lin','yscale','lin')
xlabel('ln \tau')
ylabel('\rho(ln \tau)')
figure(2)
hold on
%clf

p=linspace(0,1,1000);
plot(p,p.^(w1-1).*(1-p).^(w2-1)/beta(w1,w2))
xlabel('a')
ylabel('\rho(a)')