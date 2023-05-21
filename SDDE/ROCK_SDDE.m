function [T, Y] = ROCK_SDDE (F, G, h, t_0, Y_0, z, BackStory, rightBorder, state, damp)
%В этом методе для вычисления запаздывания в стохастической части используется 
%вычисленное методом значение в прошлом (с шумом).
    [A,B,C] = CreateRKCMethod(state, damp);
    %Определяем начальные условия
    t_zero = t_0;
    v = 0; T = [];
    rank = length(Y_0);
    Y = zeros(rank, ceil(rightBorder - t_0) / h);
    Fy = zeros(rank, ceil(rightBorder - t_0) / h);
    
    v = v + 1;
    T(v) = t_0;
    Y(:,v) = Y_0;
    Fy(:,v) = Y_0;

    k = zeros(state, rank);
    Y_1 = zeros(rank, 1);
    
    %rng('default');
    rng('shuffle');
    dW = sqrt(h)*randn(1, ceil(rightBorder-t_0)/h +1);
    
    while (t_0 < rightBorder)
         Winc = dW(v);
                
         for ii= 1:length(z)
             t_past = t_0 - z(ii);
             Delay(:,ii) = interpolation(t_past, t_zero, h, BackStory, T, Fy);
         end
%%%%%%%%%%% Высчитываем коэффициент k %%%%%%%%%%%%%%
        k(1,:) = F(t_0, Y_0, Delay);
        for i=2:state %%% Находим оставшиеся k
            y_step = Y_0;
            for j=1:rank 
                for m=1:i-1
                    y_step(j) = y_step(j) + A(i, m) * k(m, j) * h; %y_i = y_0+k*A(i,:)'*h
                end
            end
            %%%%
            for ii= 1:length(z)
               tii = t_0 + h*C(i) - z(ii);
               Delay(:,ii) = interpolation(tii, t_zero, h, BackStory, T, Fy);
            end
            k(i,:) = F( t_0 + h*C(i), y_step, Delay );
        end
        
        for j=1:rank %%% Calculed new Y
            sumY = 0;
            for i = 1:state
              sumY = sumY + h*B(i)*k(i,j);
            end
            Y_1(j) = Y_0(j) + sumY; % y_1 = y_0+K*b*h И вообще надо транспонировать k
        end
        
        Fy(:,v+1)=Y_1;
        
        for ii= 1:length(z)
            tii = t_0 + h - z(ii);
            Delay(:,ii) = interpolation(tii, t_zero, h, BackStory, T, Y);
        end
        
        Y_1 = Y_1 + G(t_0+h, Y_1, Delay)*Winc;
        
        t_0 = t_0 + h;
        Y_0 = Y_1;
        
        v = v + 1;
        T(v) = t_0;
        Y(:,v) = Y_0;
    end
end


function [A,B,C] = CreateRKCMethod (s, eta)
    A = zeros(s,s);
    B = zeros(1,s);
    T = @(s, x) Chebyshev1(s,x);
    omega0 = 1 + eta/s^2;
    omega1 = Chebyshev1(s, omega0)/Chebyshev1Diff(s, omega0);
    A(2, 1) = omega1/omega0;
    for i=3:s
            A(i, i - 1) = 2*T(i - 2, omega0)*omega1/T(i-1, omega0);
            for j=(i - 2):-1:1
                alpha = T(i - 2, omega0)/T(i-1, omega0);
                beta = T(i - 3, omega0)/T(i-1, omega0);
                A(i, j) = 2*alpha*omega0*A(i - 1, j) - beta*A(i - 2, j);
            end
    end
    
    B(s) = 2*T(s-1, omega0)*omega1/T(s, omega0);
    for j=(s - 1):-1:1
        B(j) = 2*T(s-1, omega0)*omega0*A(s, j)/ T(s, omega0) - T(s-2, omega0)*A(s - 1, j)/ T(s, omega0);
    end

    C = zeros(1,s);
    for i=1:s
        C(i)=sum(A(i,:));
    end
end