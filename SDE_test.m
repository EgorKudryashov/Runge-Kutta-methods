randn("state", 100);
lambda = 2; mu = 1; Xzero = 1;
T =1; N = 2^12; dt= T/N;

dW = sqrt(dt)*randn(1,N);
W = cumsum(dW);
Xtrue = Xzero*exp((lambda - 0.5*mu^2)*(dt:dt:T)+mu*W);
R = 16;
Dt=R*dt;
L=N/R;
Xem=zeros(1,L);
Xtemp=Xzero;
for j=1:L
    Winc = sum(dW(R*(j-1)+1:R*j));
    Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winc;
    Xem(j)=Xtemp;
end

t = linspace(0,1,length(Xtrue));
plot(t, Xtrue, 'b-');
hold on
plot([0:Dt:T], [Xzero,Xem], 'r-');
emerr = abs(Xem(end) - Xtrue(end));
hold off
