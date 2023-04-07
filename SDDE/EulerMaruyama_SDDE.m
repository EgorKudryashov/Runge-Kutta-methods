function [T, Y] = EulerMaruyama_SDDE (F, G, h, t_0, Y_0, z, BackStory, rightBorder)
    %Определяем начальные условия
    t_zero = t_0;
    v = 0; T = [];
    rank = length(Y_0);
    Y = zeros(rank, ceil(rightBorder - t_0) / h);
   
    v = v + 1;
    T(v) = t_0;
    Y(:,v) = Y_0;
    
    k = zeros(1, rank);
    Y_1 = zeros(rank, 1);
    
    rng('default');
    dW = sqrt(h)*randn(1, ceil(rightBorder-t_0)/h);
    
    while (t_0 < rightBorder)
         Winc = dW(v);
                
         for ii= 1:length(z)
             t_past = t_0 - z(ii);
             Delay(:,ii) = interpolation(t_past, t_zero, h, BackStory, T, Y);
         end
%%%%%%%%%%% Высчитываем коэффициент k %%%%%%%%%%%%%%
        k(1,:) = F(t_0, Y_0, Delay);
        
        for j=1:rank %%% Calculed new Y
            Y_1(j) = Y_0(j) + h*k(1,j); 
        end  
        Y_1 = Y_1 + G(t_0, Y_0, Delay)*Winc;
        
        t_0 = t_0 + h;
        Y_0 = Y_1;
        
        v = v + 1;
        T(v) = t_0;
        Y(:,v) = Y_0;
    end
end

