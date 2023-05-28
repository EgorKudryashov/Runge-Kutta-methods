# Runge-Kutta-methods :microscope:
Программная реализация методов Рунге-Кутты для решения задачи Коши.

![https://img.shields.io/badge/Language-MATLAB-blue](https://img.shields.io/badge/LANGUAGE-MATLAB-blue)
## Описание :scroll:
Репозиторий состоит из четырех разделов:
  1. **Корневая папка**. Здесь располагаются скрипты для графического изображения решений.
  2. **Methods**. Содержит функциональные скрипты, в которых реализованы методы решения дифференциальных уравнений.
  3. **Stability_domain**. Содержит Maple файлы, в которых строились графики областей устойчивости методов.
  4. **helpfull_function**. Дополнительные переиспользуемые функции, которые применяются в методах. 
  5. **SDDE**. Методы и скрипты посвященные решению и анализу ROCK-методов для стохастических дифференциальных уравнений с запаздыванием.
## SDDE 
### Методы :sparkle:
**EulerMaruyama_SDDE.m** - Реализация метода Эйлера-Маруямы для решения СДУЗА.

**ROCK_SDDE(2).m** - Реализация ROCK-методов для решения СДУЗА. (1) версия использует для вычисления запаздывания в шуме вычисленное значение в прошлом. (2) версия вычисляет запаздывающий шум с помощью вычисленных в прошлом значений аналогичного ДУЗА, т.е. значения функции без шума.
### Исследования :mag_right:
**ConvergeComparisonEMandROCK.m** - анализ порядка сходимости ROCK-метода. В качестве истинного решения использовалось решение методом Эйлера-Маруямы взятое с малой величиной шага.

**MeanSquareAnalisys.mlx** - анализ среднеквадратичной сходимости ROCK-методов на тестовом уравнении: TestEquationROCK.m
## Реализованные методы
### ExplicitRungeKuttaMethod.m
Реализация явных методов Рунге-Кутты с постоянным шагом.
### AutoStep_ERKMethod.m
Реализация явных методов Рунге-Кутты с автоматическим выбором шага. 
### DormanPrince.m
Реализация вложенного 7(6) этапного метода Дормана-Принца.
### ImplicitEuler.m
Реализация неявного метода Эйлера. Простейший 1-этапный неявный метод с постоянной величиной шага.
### ImplicitMidPoint.m
Реализация неявного метода по правилу средний точки. 1-этапный метод второго порядка с постоянной величиной шага.
### ImplicitHammerHollingsworth.m
Реализация неявного метода  Хаммера-Холлингсуорта. 2-этапный метод 4 порядка с постоянной величиной шага.
### DiagonalIRK.m
Реализация диагонально неявного метода. 2-этапный метод 3 порядка с постоянной величиной шага.
### RKC1.m
Реализация семейства методов Рунге-Кутты-Чебышёва для ОДУ. Можно выбрать любое количество этапов и величину демпфирования. Методы с постоянным шагом.
### RKwithDelay.m
Реализация явных методов Рунге-Кутты для решения ДУЗА. Можно выбрать метод RKC 1 порядка и задать его параметры. Функция с запаздыванием вычисляется с помощью линейной интерполяции. Методы с постоянным шагом.
### RKCforSDE.m
Реализация семейства методов Рунге-Кутты-Чебышёва для СДУ. Можно выбрать любое количество этапов и величину демпфирования. Методы с постоянным шагом. Сильный порядок p=1/2.
### EulerMaruyama.m
Классический метод Эйлера-Маруямы с постоянным шагом для СДУ.
### Milstein.m ⛔:
Недоработанная реализация метода Мильштейна для СДУ. На текущий момент работает только для автономных задач первого порядка.
