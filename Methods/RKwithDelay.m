function [X,Y] = RKwithDelay(F, x_0, y_0, z,Backstory, h, rightBorder, state, isRKC, damp)
    if (~isRKC)
        switch state %Число этапов и выбор метода
            case 1 %Метод Эйлера
                A = 0;
                B = 1;
                C = 0;
            case 2 %Метод Хойна
                A = [0 0; 1 0];
                B = [0.5 0.5];
                C = [0 1];
            case 3 %Метод Симпсона
                A = [0 0 0; 0.5 0 0; -1 2 0];
                B = [1/6 2/3 1/6];
                C = [0 0.5 1];
            case 4 %Метод Рунге-Кутты
                A = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];
                B = [1/6 1/3 1/3 1/6];
                C = [0 0.5 0.5 1];
        end
    else
        [A,B,C] = GreateRKCMethod(state, damp)
    end
    %Определяем начальные условия
    t_0 = x_0;
    v = 0; X = [];
    v = v + 1;
    X(v) = (x_0);
    Y = y_0;
    
    rank = length(y_0);
    k = zeros(state, rank);
    Y_1 = zeros(rank, 1);
    while (x_0 < rightBorder)
        %%%
         for ii= 1:length(z)
             tii = x_0 - z(ii);
             Delay(:,ii) = interpolation(tii, t_0, h, Backstory, X, Y);
         end
%         if (x_0 - z) < t_0
%             Delay = Backstory(x_0 - z);
%         else 
%             ind = bin_search(x_0-z, X);
%             theta = ((x_0-z) - X(ind)) / h;
%             Delay = Y(:,ind) * (1-theta) + Y(:,ind+1)*theta;
%         end
        %%%
%%%%%%%%%%% Высчитываем коэффициенты k %%%%%%%%%%%%%%
        k(1, :) = F(x_0, y_0, Delay);
        for i=2:state %%% Находим оставшиеся k
            y_step = y_0;
            for j=1:rank 
                for m=1:i-1
                    y_step(j) = y_step(j) + A(i, m) * k(m, j) * h; %y_i = y_0+k*A(i,:)'*h
                end
            end
            %%%%
%             for ii = 1:length(z)
%                 tii = x_0+ h*C(i) - z(ii);
%                 if (tii) < t_0
%                     Delay(:,ii) = Backstory(x_0+ h*C(i) - z);
%                 else 
%                     ind = bin_search(tii, X);
%                     theta = (tii - X(ind)) / h;
%                     Delay(:,ii) = Y(:,ind) * (1-theta) + Y(:,ind+1)*theta;
%                 end
%             end
           %%%%
            for ii= 1:length(z)
               tii = x_0 + h*C(i) - z(ii);
               Delay(:,ii) = interpolation(tii, t_0, h, Backstory, X, Y);
            end
            k(i,:) = F( x_0 + h*C(i), y_step, Delay );
        end
        for j=1:rank %%% Calculed new Y
            sum = 0;
            for i = 1:state
              sum = sum + h*B(i)*k(i,j);
            end
            Y_1(j) = y_0(j) + sum; % y_1 = y_0+K*b*h И вообще надо транспонировать k
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        x_0 = x_0 + h;
        y_0 = Y_1;
        
        v = v + 1;
        X(v) = (x_0);
        Y = [Y, y_0];
    end
   
end

function [A,B,C] = GreateRKCMethod (s, eta)
    A = zeros(s,s);
    T = @(s, x) Chebyshev1(s,x);
    omega0 = 1 + eta/s^2;
     omega1 = Chebyshev1(s, omega0)/Chebyshev1Diff(s, omega0);
%     if (eta==0)
%         omega1 = 1/s^2;
%     else
%         omega1 = T(s, omega0) / (s * (T(s-1, omega0) -  omega0 * T(s, omega0))) * (1 - omega0^2);
%     end
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
    sum(B)
%     if (sum(B)<1)
%         B(1)= B(1)+1-sum(B);
%     end
    C = zeros(1,s);
    for i=1:s
        C(i)=sum(A(i,:));
    end
end
