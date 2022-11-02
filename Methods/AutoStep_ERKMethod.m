function [X,Y] = AutoStep_ERKMethod(F, state, x_0, y_0, h, rightBorder, tol)
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
    
    %Определяем начальные условия
    p = state; %порядок сходимости
    v = 0; X = []; % rank= length(F); k = zeros(state, rank); Y_1 = zeros(1,rank);
    v = v + 1;
    X(v) = (x_0);
    Y = y_0;
    
    while (x_0 < rightBorder)
        Y_1 = Step (F, A,B,C, h, x_0, y_0);
    %%%Высчитываем интервал половинчатым шагом%%%%%%%%
        h_halfStep = h/2;
        y_halfStep = y_0; x_halfStep = x_0; 

        Y_1_halfStep = Step(F,A,B,C, h_halfStep, x_halfStep, y_halfStep);
        Y_polovinka = Y_1_halfStep;
        x_halfStep = x_halfStep + h_halfStep;
        Y_1_halfStep = Step(F,A,B,C, h_halfStep, x_halfStep, Y_polovinka);

    %%%%Находим погрешность
        r_n = norm(abs(Y_1 - Y_1_halfStep)) / ((2^(p)-1));

        %%%%Делаем выбор%%%%%%%%%%%%%%%%%%%%%%%%%%
        v = v + 1;
        %1 условие, плохо, пересчитываем шаг 
        while ( r_n > (tol*2^p) )
            h = h_halfStep;
            Y_1 = Y_polovinka;

            h_halfStep = h/2;
            y_halfStep = y_0;
            x_halfStep = x_0;

            %%Высчитываем интервал половинчатым шагом%%%%%%%%        
            Y_1_halfStep = Step(F,A,B,C, h_halfStep, x_halfStep, y_halfStep);
            Y_polovinka = Y_1_halfStep;
            x_halfStep = x_halfStep + h_halfStep;
            Y_1_halfStep = Step(F,A,B,C, h_halfStep, x_halfStep, Y_polovinka);

            %%Находим погрешность
            r_n = norm(abs(Y_1 - Y_1_halfStep)) / ((2^(p)-1));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %2 условие, норм, но принимаем меньший шаг
        if (tol<r_n) && (r_n<=tol*2^p) 
            x_0 = x_0 + h;
            y_0 = Y_1_halfStep;

            h = h_halfStep;
        end

        %3 условие, супер, шаг оставляем
        if (tol*(1/2^(p+1))<=r_n) && (r_n<=tol)
            x_0 = x_0 + h;
            y_0 = Y_1;
        end

        %4 условие,очень точно, увеличиваем шаг
        if (r_n<tol*(1/(2^(p+1))))    
            x_0 = x_0 + h;
            y_0 = Y_1;

            h = h * 2;
        end

        X(v) = (x_0);
        Y = [Y, y_0];
        
    end
    
end

