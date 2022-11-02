%%
% Описание хим. реакций
F = @(x,y)[ -0.04*y(1) + 10^4 * y(2) * y(3); 
             0.04*y(1) - 10^4 * y(2) * y(3) - 3*10^7* y(2)*y(2);
             3*10^7* y(2)*y(2)];
x_0 = 0;
y_0 = [1; 0; 0];

RightBorder = 4;
h = 1/5;
tol = 10^(-8);
% xspan = [0, 4];
% [x,f]= ode45(F, xspan, y_0);

[X,Y] = ImplicitHammerHollingsworth(F, x_0,y_0, h, RightBorder);
[Xd,Yd] = DormanPrince(F, x_0,y_0,h/100,RightBorder, tol);

% figure
% hold on
% grid on
% plot (X,Y(1,:),'-o',x,f(:,1));
% plot (X,Y(2,:),'-o',x,f(:,2));
% plot (X,Y(3,:),'-o',x,f(:,3));

figure
hold on
grid on
plot (X,Y(1,:),'-o',Xd,Yd(1,:));
plot (X,Y(2,:),'-o',Xd,Yd(2,:));
plot (X,Y(3,:),'-o',Xd,Yd(3,:));
title ('tol = 10^{-6}, Число шагов = 1986 для метода Дорманда-Принса')
legend ('Неявный метод с 2-стадиями', 'Метод Дорманда-Принса')
xlabel('x')
ylabel('y')
%%
%Простой пример жесткой задачи
F = @(x,y) -50*(y - cos(x));
x_0 = 0;
y_0 = 0;
h = 1/25;
RightBorder = 2;

[X,Y] = ExplicitRungeKuttaMethod(F, 1, x_0, y_0, h, RightBorder);
[Xim, Yim] = ImplicitEuler (F, x_0, y_0, 0.5, RightBorder);
[x,y] = DormanPrince(F, x_0, y_0, h, RightBorder, 10^-6);
figure
hold on
grid on
plot (X,Y,Xim,Yim,'-o');
plot(x,y, 'r');
title ('h = 1/25 для явного метода Эйлера, h = 1/2 для неявного')
legend ('явный метод Эйлера', 'неявный метод Эйлера', 'вложенный метод Дорманда-Принса')

%%
F = @(x,y) [y(2);               %уравнение Ван-Дер Поля
           20*((1-y(1)^2) * y(2)) - y(1)];
x_0 = 0;
y_0 = [2; 0];
RightBorder = 16;

%[X1,Y1]= DiagonalIRK (F, x_0, y_0, 1/20, RightBorder);
[X1,Y1]= ImplicitEuler (F, x_0, y_0, 1, RightBorder);
[X2,Y2]=ImplicitHammerHollingsworth(F, x_0,y_0, 1, RightBorder);
[x,y] = ODIRK (F, x_0,y_0, 1/100, RightBorder, 10^-6);
%[x,y] = ExplicitRungeKuttaMethod(F, 4, x_0, y_0, 1/20, RightBorder);

figure

hold on
grid on
plot (Y1(1,:),Y1(2,:),'-*');
plot (Y2(1,:),Y2(2,:),'-o');
plot (y(1,:), y(2,:),'green')
xlim([-3, 3]);
ylim([-RightBorder, RightBorder]);
ylabel('y_2')
xlabel ('y_1')
title ('h = 1/20 для явного метода, h = 1/2 для неявных')
%legend ('Диагонально неявный метод 3-го порядка', 'Неявный метод Хаммера-Холлингсуорта 4-го порядка', 'Явный метод 4-го порядка')
legend ('Неявный метод Эйлера', 'Неявный метод Хаммера-Холлингсуорта 4-го порядка', 'Точное решение')

figure
tiledlayout(2,1);
nexttile;
hold on
grid on

plot (X1(1,:),Y1(1,:),'-*');
plot (X2(1,:),Y2(1,:),'-o');
plot (x(1,:), y(1,:), 'green')
xlim([0, RightBorder]);
ylim([-3, 3]);
xlabel('x')
ylabel ('y_1')
title ('h = 1 для неявных методов')
%title ('h = 1/20 для всех методов')
%legend ('Диагонально неявный метод 3-го порядка', 'Неявный метод Хаммера-Холлингсуорта 4-го порядка', 'Явный метод 4-го порядка')
legend ('Неявный метод Эйлера', 'Неявный метод Хаммера-Холлингсуорта 4-го порядка', 'Точное решение')

nexttile;
hold on
grid on

plot (X1(1,:),Y1(2,:),'-*');
plot (X2(1,:),Y2(2,:),'-o');
plot (x(1,:), y(2,:), 'green')
xlim([0, RightBorder]);
ylim([-40, 40]);
xlabel('x')
ylabel ('y_2')
title ('h = 1 для неявных методов')
%title ('h = 1/20 для всех методов')
%legend ('Диагонально неявный метод 3-го порядка', 'Неявный метод Хаммера-Холлингсуорта 4-го порядка', 'Явный метод 4-го порядка')
legend ('Неявный метод Эйлера', 'Неявный метод Хаммера-Холлингсуорта 4-го порядка', 'Точное решение')