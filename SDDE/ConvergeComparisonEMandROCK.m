clear all;

rng('shuffle');

%%% Task description
lambda = -2; mu = -1;
F = @(t,y,z) lambda*y + 0.5*z;  G = @(t,y,z) mu*y + 0.2*z;
history = @(t,y) 1;
delay = 1/4;
T_0 = 0;
Y_0 = 1;

rightBorder = 1; N = 2^18; dt= rightBorder/N;

setH = 5;
line = zeros(1,setH);
for p=1:setH
   R(p) = 2^(setH+10-p);
   line(p)= R(p)*dt;
   intervals(p) = ceil(N/R(p)) + 1;
end

% Description ROCK method
state = 4;
damp = 1;
[A,B,C] = CreateRKCMethod(state, damp);
%A = [0 0; 0.2785 0];
%B = [0.5 0.5];
%C = [0 0.2785];

    
numExperiments = 1000;
WEAK = zeros(1,setH);
STRONG = zeros(1,setH);

ALLTrue = zeros(numExperiments, N+1);
ALLROCK = zeros(numExperiments, setH, max(intervals));
for exp=1:numExperiments

%%% 'True' solution
    dW = sqrt(dt)*randn(1,N);

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
    
    ALLTrue(exp, :) = Ytrue;

    for p = 1:setH
        
        Dt=R(p)*dt;
        L=ceil(N/R(p));
        
%%%%%%%%%% ROCK solution
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
        
        %%% Для ROCK_SDDE1 - закомментировать и заменить на Yrock в вызове
        %%% интерполяции для delay в шуме
        Fy = zeros(rank, ceil(L));
        Fy(:,v) = y_0;

        while (t_0 < rightBorder)
            Winc = 0;
            for l = R(p)*(v-1)+1: R(p)*v
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
                 Y_1(j) = y_0(j) + sumY;
             end
            %%%
            Fy(:,v+1)=Y_1;
            %%%
             for ii= 1:length(delay)
                 tii = t_0 + Dt - delay(ii);
                 Z(:,ii) = interpolation(tii, T_0, Dt, history, Trock, Fy);
             end

            Y_1 = Y_1 + G(t_0+Dt, Y_1, Z)*Winc;

            v = v+1;
            t_0 = t_0 + Dt;
            y_0 = Y_1;

            Trock(v) = t_0;
            Yrock(:,v) = y_0;
        end

        for i2=1:intervals(p)
            ALLROCK(exp, p, i2) = Yrock(i2);
        end
    end
end

for p=1:setH
    for i2=1:intervals(p)
        tmpWeak(i2) = abs( mean(ALLTrue(:, 1+(i2-1)*R(p))) - mean(ALLROCK(:,p,i2)) ); 
    end
    WEAK(p) = max(tmpWeak);
end

for p=1:setH
    for i2=1:intervals(p)
       tmpStrong(i2) = mean(abs( ALLTrue(:, 1+(i2-1)*R(p)) - ALLROCK(:,p,i2) ));
    end
    STRONG(p) = max(tmpStrong);
end
   
%%% Plot building
figure
plot(Ttrue, Ytrue, 'b-');
hold on
plot(Trock, Yrock, 'r-');
legend("Точное решение","ROCK")
hold off


figure
subplot(1,2,1)
loglog(line, WEAK, '*-');
hold on
loglog(line, line / line(1) * WEAK(1), '--');
hold off
title("Слабая сходимость");
xlabel("h");
ylabel("max |E[ y_n ] - E[ y( t_n ) ]|");
ax = gca;
ax.FontSize = 14;
legend("Слабая ошибка для ROCK","Наклон 1/2")

subplot(1,2,2)
loglog(line, STRONG, '*-');
hold on
loglog(line, line.^(1/2) / sqrt(line(1)) * STRONG(1), '--');
hold off
title("Сильная сходимость");
xlabel("h");
ylabel("max E[| y_n - y( t_n ) |]");
ax = gca;
ax.FontSize = 14;
legend("Сильная ошибка для ROCK","Наклон 1")
