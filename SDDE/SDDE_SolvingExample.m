clear all;

rng('shuffle');

%%% Task description
lambda = -2; mu = -1;
F = @(t,y,z) lambda*y+z(1,1);  G = @(t,y,z) mu*y+z(1,2);
history = @(t,y) 1;
delay = [1/2, 1/4];
T_0 = 0;
Y_0 = 1;

rightBorder = 6; N = 2^16; dt= rightBorder/N;

setH = 10;
line = zeros(1,setH);
WEAK = zeros(1,setH);
STRONG = zeros(1,setH);

for p = 1:setH
    
dW = sqrt(dt)*randn(1,N);

%%% Parameters for working with convergences
numExperiments = 100;
YmeanTrue = zeros(numExperiments, N+1);

R = 2^(setH+2-p);
Dt=R*dt;
L=ceil(N/R);

line(p) = Dt;

YmeanSolution = zeros(numExperiments, L+1);
AverageSolution = zeros(numExperiments, L+1);


% Description ROCK method
state = 2; damp = 0.1;
A = [0 0; 0.2785 0];
B = [0.5 0.5];
C = [0 0.2785];


%%% 'True' solution
for exp=1:numExperiments
    
    v = 0; Ttrue = [];
    t_0 = T_0;
    y_0 = Y_0;
    rank = length(y_0);
    Ytrue = zeros(rank, N);
   
    v = v + 1;
    Ttrue(v) = t_0;
    Ytrue(:,v) = y_0;
    k = zeros(1, rank);
    Y_1 = zeros(rank, 1);
  
    while (t_0 < rightBorder)
        Winc = dW(v);             
        for ii= 1:length(delay)
            t_past = t_0 - delay(ii);
            Z(:,ii) = interpolation(t_past, T_0, dt, history, Ttrue, Ytrue);
        end
        k(1,:) = F(t_0, y_0, Z);
        
        for j=1:rank %%% Calculed new Y
            Y_1(j) = y_0(j) + dt*k(1,j); 
        end  
        Y_1 = Y_1 + G(t_0, y_0, Z)*Winc;
        
        t_0 = t_0 + dt;
        y_0 = Y_1;
        
        v = v + 1;
        Ttrue(v) = t_0;
        Ytrue(:,v) = y_0;
    end

    YmeanTrue(exp,:) = Ytrue;

%%% ROCK solution
    v = 0; Trock = [];
    t_0 = T_0;
    y_0 = Y_0;
    rank = length(y_0);
    Yrock = zeros(rank, ceil(L));
   
    v = v + 1;
    Trock(v) = t_0;
    Yrock(:,v) = y_0;
    k = zeros(1, rank);
    Y_1 = zeros(rank, 1);
    
while (t_0 < rightBorder)
    Winc = 0;
    for l = R*(v-1)+1: R*v
        if l < length(dW)
            Winc = Winc + dW(l);
        end
    end
    
    for ii= 1:length(delay)
        t_past = t_0 - delay(ii);
        Z(:,ii) = interpolation(t_past, T_0, Dt, history, Trock, Yrock);
    end
        k(1,:) = F(t_0, y_0, Z);
        for i=2:state %%% Находим оставшиеся k
            y_step = y_0;
            for j=1:rank 
                for m=1:i-1
                    y_step(j) = y_step(j) + A(i, m) * k(m, j) * Dt; %y_i = y_0+k*A(i,:)'*h
                end
            end
            %%%%
            for ii= 1:length(delay)
               tii = t_0 + Dt*C(i) - delay(ii);
               Z(:,ii) = interpolation(tii, T_0, Dt, history, Trock, Yrock);
            end
            k(i,:) = F( t_0 + Dt*C(i), y_step, Z );
        end
        
        for j=1:rank %%% Calculed new Y
            sumY = 0;
            for i = 1:state
              sumY = sumY + Dt*B(i)*k(i,j);
            end
            Y_1(j) = y_0(j) + sumY; % y_1 = y_0+K*b*h И вообще надо транспонировать k
        end
        
        for ii= 1:length(delay)
            tii = t_0 + Dt - delay(ii);
            Z(:,ii) = interpolation(tii, T_0, Dt, history, Trock, Yrock);
        end
        
        Y_1 = Y_1 + G(t_0+Dt, Y_1, Z)*Winc;
        
        v = v+1;
        t_0 = t_0 + Dt;
        y_0 = Y_1;
        
        Trock(v) = t_0;
        Yrock(:,v) = y_0;
end

    YmeanSolution(exp,:) = Yrock;
    
    for expStrong=1:L+1
        AverageSolution(exp,expStrong) = abs(Ytrue(1+(expStrong-1)*R) - Yrock(expStrong));
    end
end

% Means result
meanYtrue = zeros(1, N+1);
meanYsolving = zeros(1, L+1);
weakConv = zeros (1, L+1);
strongConv = zeros(1, L+1);

for exp=1:N+1
    meanYtrue(exp) = mean(YmeanTrue(:,exp).^2);
end

for exp=1:L+1
    meanYsolving(exp) = mean(YmeanSolution(:,exp).^2);
end

for exp=1:L+1
   weakConv(exp) = abs(meanYtrue(1+(exp-1)*R) - meanYsolving(exp));
end

for exp=1:L+1
   strongConv(exp) = mean(AverageSolution(:,exp).^2);
end

WEAK(p) = weakConv(L-1);
STRONG(p) = strongConv(L-1);
end


%%% Plot building
figure
plot(Ttrue, Ytrue, 'b-');
hold on
plot(Trock, Yrock, 'r-');
legend("Точное решение","ROCK")
hold off

figure
semilogy(Ttrue, meanYtrue, 'b-');
hold on
semilogy(Trock, meanYsolving, 'r-');
hold off

figure
subplot(2,1,1)
loglog(line, WEAK, '*-');
title("Weak convergence");

subplot(2,1,2)
loglog(line, STRONG, '*-');
title("Strong convergence");
