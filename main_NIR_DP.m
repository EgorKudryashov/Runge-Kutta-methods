F = @(x,y) [ exp(y(2));
             x*sin(y(1)/10)];
x_0 = 0;
y_0 = [0;0];
RightBorder = 4.5;
h = 0.05;
%%
[X,Y] = ExplicitRungeKuttaMethod (F, 4, x_0, y_0, h, RightBorder);
[x,f] = DormanPrince(F, x_0,y_0, h, RightBorder, 10^-6);

figure
tiledlayout(1,1);

nexttile;
hold on
plot (X, Y(1,:),'-o',x,f,'-*')
plot (X, Y(2,:),'-o',x,f,'-*')
%title("Решение системы: y'_1 = exp(y_2) and y'_2 = x*sin(y_1/10) , y_1(0)= 0, y_2(0)=0")
title("Величина шага для метода Рунге-Кутты h = 0.05")
legend ({'Классический метод Рунге-Кутты', 'Вложенный метод'},...
    'Location','northwest')
grid on
xlabel('x')
ylabel('y')
%%
%[X,Y] = DormanPrince(F, x_0,y_0, h, RightBorder, 10^-4);
[X2,Y2] = DormanPrince(F, x_0,y_0, h, RightBorder, 10^-6);
[X3,Y3] = DormanPrince(F, x_0,y_0, h, RightBorder, 10^-10);
[x,f] = DormanPrince(F, x_0,y_0, h, RightBorder, 10^-12);

figure
tiledlayout(1,2);

% nexttile;
% hold on
% plot (X, Y(1,:),'-o',x,f)
% plot (X, Y(2,:),'-o',x,f)
% title('tol = 10^{-4}, Число шагов = 3')
% legend ({'Метод Дормана-Принца', 'Точное решение'},...
%     'Location','northwest')
% grid on
% xlabel('x')
% ylabel('y')

nexttile;
hold on
plot (X2, Y2(1,:),'-o',x,f)
plot (X2, Y2(2,:),'-o',x,f)
title('tol = 10^{-6}, Число шагов = 27')
legend ({'Метод Дорманда-Принса', 'Точное решение'},...
    'Location','northwest')
grid on
xlabel('x')
ylabel('y')

nexttile;
hold on
plot (X3, Y3(1,:),'-o',x,f)
plot (X3, Y3(2,:),'-o',x,f)
title('tol = 10^{-10}, Число шагов = 194')
legend ({'Метод Дорманда-Принса', 'Точное решение'},...
    'Location','northwest')
grid on
xlabel('x')
ylabel('y')

sgtitle("Решение системы: y'_1 = exp(y_2) and y'_2 = x*sin(y_1/10) , y_1(0)= 0, y_2(0)=0")
