k0 = 1;
close all
fs = 18;
x0 = [ 1 2 2.5 5]*k0;
c = 20;
x = [2.718281828	7.389056099	12.18249396	148.4131591]*k0;
xx0 = [1:.01:5]*k0;
N = numel(xx);


x = 1./x0;
xx = 1./xx0;
X = [ones(4,1), x'];
XX = [ones(N,1), xx'];
y = [84.19	49.44	41.56	25.48]';
b = (X'*X)^-1*X'*y;
c = b(1);
y2 = log(y-c);
b2 = (X'*X)^-1*X'*y2;

figure
hold on
plot(x0,y,'o')
plot(xx0,XX*b)
plot(xx0,exp(XX*b2)+c)
legend({'Simulated Data','1/x Fit','Exponential Fit'},'FontSize',fs)

xlabel('k (k_0)','FontSize',fs)
ylabel('Median Waiting Time (a.u.)','FontSize',fs)