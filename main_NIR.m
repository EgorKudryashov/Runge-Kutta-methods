F = @(x,y) cos(x) - y;
Realf = @(x) 1/2*(-exp(-x) + sin(x) + cos(x));
x = linspace(0, 4);
f = Realf(x);
x_0 = 0;
y_0 = 0;
RightBorder = 4;
h = 1/2;

%% Метод Эйлера
[X,Y] = ExplicitRungeKuttaMethod (F, 1, x_0, y_0, h, RightBorder);
[X2,Y2] = ExplicitRungeKuttaMethod (F, 1, x_0, y_0, h/5, RightBorder);
[X3,Y3] = ExplicitRungeKuttaMethod (F, 1, x_0, y_0, h/20, RightBorder);

figure
tiledlayout(3,1);

nexttile;
plot (X, Y(1,:),'-o', x, f)
title('h=0.5, Число шагов = 8')
legend ('Метод Эйлера', 'Точное решение')
grid on
xlabel('x')
ylabel('y')

nexttile
plot (X2, Y2(1,:),'-o', x, f)
title('h=0.1, Число шагов = 40')
legend ('Метод Эйлера', 'Точное решение')
grid on
xlabel('x')
ylabel('y')
hold on
grid on

nexttile
plot (X3, Y3(1,:),'-o', x, f)
title('h=0.025, Число шагов = 160')
legend ('Метод Эйлера', 'Точное решение')
grid on
xlabel('x')
ylabel('y')
hold on
grid on

%% Методы Хойна и Симпсона
[X,Y] = ExplicitRungeKuttaMethod (F, 2, x_0, y_0, h, RightBorder);
[X2,Y2] = ExplicitRungeKuttaMethod (F, 2, x_0, y_0, h/5, RightBorder);

[Xs,Ys] = ExplicitRungeKuttaMethod (F, 3, x_0, y_0, h, RightBorder);
[X2s,Y2s] = ExplicitRungeKuttaMethod (F, 3, x_0, y_0, h/5, RightBorder);

figure
tiledlayout(2,2);

nexttile;
plot (X, Y(1,:),'-o', x, f)
title('h=0.5, Число шагов = 8')
legend ('2-этапный метод', 'Точное решение')
grid on
xlabel('x')
ylabel('y')

nexttile
plot (X2, Y2(1,:),'-o', x, f)
title('h=0.1, Число шагов = 40')
legend ('2-этапный метод', 'Точное решение')
grid on
xlabel('x')
ylabel('y')
hold on
grid on

nexttile;
plot (Xs, Ys(1,:),'-o', x, f)
title('h=0.5, Число шагов = 8')
legend ('3-этапный метод', 'Точное решение')
grid on
xlabel('x')
ylabel('y')

nexttile
plot (X2s, Y2s(1,:),'-o', x, f)
title('h=0.1, Число шагов = 40')
legend ('3-этапный метод', 'Точное решение')
grid on
xlabel('x')
ylabel('y')
hold on
grid on

%% Классический метод Рунге-Кутты

F = @(x,y) cos(x) - y;
Realf = @(x) 1/2*(-exp(-x) + sin(x) + cos(x));
x = linspace(0, 4);
f = Realf(x);
x_0 = 0;
y_0 = 0;
RightBorder = 4;
h = 1;

[X,Y] = ExplicitRungeKuttaMethod (F, 4, x_0, y_0, h, RightBorder);

figure
tiledlayout(2,1);

nexttile;
plot (X, Y(1,:),'-o', x, f)
title('h= 1, Число шагов = 4')
legend ('Метод Рунге-Кутты', 'Точное решение')
grid on
xlabel('x')
ylabel('y')


F = @(x,y) [ x*y(2)^(1/4)*y(4);
             x*exp(y(3)+2)*y(4)^(2);
             sin(y(4));
             -1*x*log(y(1));
             ];
x_0 = 0;
y_0 = [1; 1; -2; 1];
RightBorder=4;
h = 0.05;

[x,f] = DormanPrince(F,x_0,y_0, h, RightBorder, 10^-10);
[X,Y] = ExplicitRungeKuttaMethod (F, 4, x_0, y_0, h, RightBorder);

nexttile;
hold on
plot (X, Y(1,:),'-o', x, f(1,:))
plot (X, Y(2,:),'-o', x, f(2,:))
plot (X, Y(3,:),'-o', x, f(3,:))
plot (X, Y(4,:),'-o', x, f(4,:))
title('h= 0.05, Число шагов = 80')
legend ('Метод Рунге-Кутты', 'Точное решение')

grid on
xlabel('x')
ylabel('y')

%% График сравнения всех явных методов
[X,Y] = ExplicitRungeKuttaMethod (F, 1, x_0, y_0, h/2, RightBorder);
[X2,Y2] = ExplicitRungeKuttaMethod (F, 2, x_0, y_0, h, RightBorder);
[X3,Y3] = ExplicitRungeKuttaMethod (F, 3, x_0, y_0, h, RightBorder);
[X4,Y4] = ExplicitRungeKuttaMethod (F, 4, x_0, y_0, h*2, RightBorder);

figure

tiledlayout(2,2);


nexttile;
plot (X, Y(1,:),'-o', x, f)
title('h = 0.25, Число шагов = 16')
legend ('Метод 1-го порядка', 'Точное решение')
grid on
xlabel('x')
ylabel('y')

nexttile
plot (X2, Y2(1,:),'-o', x, f)
title('h = 0.5, Число шагов = 8')
legend ('Метод 2-го порядка', 'Точное решение')
grid on
xlabel('x')
ylabel('y')
hold on
grid on

nexttile
plot (X3, Y3(1,:),'-o', x, f)
title('h = 0.5, Число шагов = 8')
legend ('Метод 3-го порядка', 'Точное решение')
grid on
xlabel('x')
ylabel('y')
hold on
grid on

nexttile
plot (X4, Y4(1,:),'-o', x, f)
title('h = 1, Число шагов = 8')
legend ('Метод 4-го порядка', 'Точное решение')
grid on
xlabel('x')
ylabel('y')
hold on
grid on

sgp = sgtitle (" Решение задачи: y' = cos(x) - y,  y(0) = 0 ", 'Color', 'black');
sqp.FontSize = 24;