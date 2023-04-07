function [T, Y] = RKCforSDE(F, G, h, t_0, Y_0, rightBorder, state, damp)

    %[A,B,C] = GreateRKCMethod(state, damp);
    A = [0 0; 0.2785 0];
    B = [0.5 0.5];
    C = [0 0.2785];
    T = zeros(1, ceil(rightBorder - t_0) / h);
    rank = length(Y_0);
    Y = zeros(rank, ceil(rightBorder - t_0) / h);
    v = 1;
    T(v) = t_0;   Y(:,v) = Y_0;
    Ytemp = Y_0;

    %rng('default');
    rng('shuffle');
    m = 1;
    dW = sqrt(h/m)*randn(1, ceil(m*(rightBorder-t_0)/h));
    
    k = zeros(state, rank);
    while t_0 < rightBorder
        Winc = sum(dW((m*(v-1)+1):(m*v)));
        
        k(1, :) = F(t_0, Ytemp);
        for i=2:state %%% Находим оставшиеся k
            y_step = Ytemp;
            for j=1:rank 
                for l=1:i-1
                    y_step(j) = y_step(j) + A(i, l) * k(l, j) * h; %y_i = y_0+k*A(i,:)'*h
                end
            end
            k(i,:) = F( t_0 + h*C(i), y_step);
         end
        for j=1:rank %%% Calculed new Y
             sumY = 0;
             for i = 1:state
                sumY = sumY + h*B(i)*k(i,j);
             end
             Ytemp(j) = Ytemp(j) + sumY; % y_1 = y_0+K*b*h И вообще надо транспонировать k
        end
        
        Ytemp = Ytemp + G(t_0+h, Ytemp)*Winc;
    
        t_0 = t_0 + h;
        v = v+1;
        T(v) = t_0;
        Y(:,v) = Ytemp;
    end
end

function [A,B,C] = GreateRKCMethod (s, eta)
    A = zeros(s,s);
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
