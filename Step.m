function [Y_1] = Step(F, A, B, C, h, x_0, y_0 )
    rank = length(y_0);
    state = length(B);
    k = zeros(state, rank);
 %%%% Высчитываем коэффициенты k %%%%%%%%%%%%%
        for i=1:rank  %%% Находим k_1
            k(1, i) = F{i}( x_0, y_0 );
        end
        for i=2:state %%% Находим оставшиеся k
            y_step = y_0;
            for j=1:rank 
                for m=1:i-1
                    y_step(j) = y_step(j) + A(i, m) * k(m, j) * h;
                end
            end
            for j=1:rank
                k(i,j) = F{j}( x_0 + h*C(i), y_step );
            end
        end
        for j=1:rank %%% Считаем новые Y
            sum = 0;
            for i = 1:state
              sum = sum + h*B(i)*k(i,j);
            end
            Y_1(j) = y_0(j) + sum;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end

