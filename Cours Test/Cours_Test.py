#метод Рунге—Кутта—Фельберга с автоматическим выбором шага 
from numpy import*
import matplotlib.pyplot as plt
import numpy as np
import time




def RungeKuttaFehlborga(f, to, yo, zo, tEnd, h, eps):
		def increment(f, t, y, h, flag = 0):# поиск приближённого решения методом Рунге—Кутта—Фельберга.
			k1=h*f(t,y)
			k2=h*f(t+(1/4)*h,y+(1/4)*k1)
			k3=h*f(t+(3/8)*h,y+(3/32)*k1+(9/32)*k2)
			k4=h*f(t+(12/13)*h,y+(1932/2197)*k1-(7200/2197)*k2+(7296/2197)*k3)
			k5=h*f(t+h,y+(439/216)*k1-8*k2+(3680/513)*k3 -(845/4104)*k4)
			k6=h*f(t+(1/2)*h,y-(8/27)*k1+2*k2-(3544/2565)*k3 +(1859/4104)*k4-(11/40)*k5)
			if flag == 4:
				return (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - k5/5 #4 порядка 
			else:
				return (16/135)*k1+(6656/12825)*k3+(28561/56430)*k4-(9/50)*k5+(2/55)*k6 #5 порядка
		t = []#подготовка пустого списка t
		y = []#подготовка пустого списка y
		z = []#подготовка пустого списка z

		t.append(to)#внесение в список t начального значения to
		y.append(yo)#внесение в список y начального значения yo
		z.append(zo)#внесение в список y начального значения zo
		while to < tEnd:#внесение результатов расчёта в массивы t,y
			
			if (abs(yo - zo) > eps):
				h = h/2
			elif (abs(yo - zo) < eps/32):
				h = h*2
			elif (eps/32 < abs(yo - zo) < eps):
				h = h

			zo = yo + increment(f, to, yo, h, 4)
			yo = yo + increment(f, to, yo, h) # расчёт значения в точке t0,y0 для задачи Коши
			to = to + h # приращение времени
			t.append(to) # заполнение массива t
			y.append(yo) # заполнение массива y   
			z.append(zo) # заполнение массива z
		return array(t), array(y), array(z)


def f(t,y):
	return math.cos(t)*t**2 + math.sin(t)*2*t

to = 0
tEnd = 50
yo = 0
h = 1
eps = 0.000001

t, y, z = RungeKuttaFehlborga(f, to, yo, yo, tEnd, h, eps)
for i in range(t.size):
	print(str(t[i]))
	print(str(y[i]))
	print(str(z[i]))


plt.title('Результат численного решение ОДУ 1-го порядка')
plt.plot(t, z, '-', color = 'red')	
plt.plot(t, y, '--', color = 'blue')
plt.grid()
plt.show()


f = open("values.txt", "w")
for i in range(t.size):
	f.write("t(" + str(i) + "): " + str(t[i]) + "      " +"y(" + str(i) + "): " +  str(y[i]) + '\n')
f.close()