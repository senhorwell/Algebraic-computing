import math

# Put here the function that has to be calculated
 
f  = lambda x, y: 2 - math.exp(-4 * x) - 2 * y

# Here is the information of the exercise

x0 = 0
y0 = 1
h  = 0.5
xn = 1

# Put here the order of the method of Runge Kutta
ordem_do_metodo = 4

# Ready! Now run "python ruge.py" on the terminal and be happy :)

# ================ FUNCOES ================

def calcK(func, x, y, h, k):
	return func(x + h, y + h * k)

def ruge2O(funcao, x0: float, y0: float, h: float, x_final: float):
	n = int((x_final - x0) / h)

	xn = x0
	yn = y0

	print("Runge-Kutta Method 2nd Order:")

	result = ["xn", "yn", "k1", "k2"]
	print(' '.join(i.rjust(10) for i in result))
	
	for _ in range(n + 1):
		
		k1 = calcK(funcao, xn, yn, 0, 0)
		k2 = calcK(funcao, xn, yn, h, k1)

		result = [xn, yn, k1, k2]
		print(' '.join('{:3.6f}'.format(i).rjust(10) for i in result))

		yn = yn + h * (0.5 * k1 + 0.5 * k2)
		xn += h

def ruge4O(funcao, x0: float, y0: float, h: float, x_final: float):
	n = int((x_final - x0) / h)

	xn = x0
	yn = y0

	print("Runge-Kutta Method 4th Order:")

	result = ["xn", "yn", "k1", "k2", "k3", "k4"]
	print(' '.join(i.rjust(10) for i in result))

	for _ in range(n + 1):
		
		k1 = calcK(funcao, xn, yn, 0,     0)
		k2 = calcK(funcao, xn, yn, h / 2, k1)
		k3 = calcK(funcao, xn, yn, h / 2, k2)
		k4 = calcK(funcao, xn, yn, h,     k3)

		result = [xn, yn, k1, k2, k3, k4]

		print(' '.join("{:3.6f}".format(i).rjust(10) for i in result))

		yn = yn + h * (
			(1 / 6) * k1 +
			(1 / 3) * k2 +
			(1 / 3) * k3 +
			(1 / 6) * k4)
		xn += h

if __name__ == "__main__":

	if ordem_do_metodo == 2:
		ruge2O(f, x0, y0, h, xn)
	elif ordem_do_metodo == 4:
			ruge4O(f, x0, y0, h, xn)