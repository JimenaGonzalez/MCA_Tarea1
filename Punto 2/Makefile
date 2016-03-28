all:12R.out plot12R.o 12L.out plot12L.o 13R.out plot13R.o 13L.out plot13L.o 9R.out plot9R.o 9L.out plot9L.o 10R.out plot10R.o 10L.out plot10L.o 11R.out plot11R.o  plot17R.o 11L.out plot11L.o plot17L.o
	




12R.out:  RK4.c
	gcc RK4.c -lm

plot12R.o:plot.py a.out
	python plot.py 'condiciones12.txt' -3 3 -3 3 RK4_12

12L.out:  simplect.c
	gcc simplect.c -lm

plot12L.o:plot.py a.out
	python plot.py 'condiciones12.txt' -3 3 -3 3 Simplectico_12

13R.out:  RK4.c
	gcc RK4.c -lm

plot13R.o:plot.py a.out
	python plot.py 'condiciones12.txt' -1 1 -1 1 RK4_13

13L.out:  simplect.c
	gcc simplect.c -lm

plot13L.o:plot.py a.out
	python plot.py 'condiciones12.txt' -1 1 -1 1. Simplectico_13



9R.out:  RK4.c
	gcc RK4.c -lm

plot9R.o:plot.py a.out
	python plot.py 'condiciones9.txt' -3 3 -2.5 2.5 RK4_9

9L.out:  simplect.c
	gcc simplect.c -lm

plot9L.o:plot.py a.out
	python plot.py 'condiciones9.txt' -3 3 -2.5 2.5 Simplectico_9


10R.out:  RK4.c
	gcc RK4.c -lm

plot10R.o:plot.py a.out
	python plot.py 'condiciones9.txt' -0.2 0.2 -0.2 0.2 RK4_10

10L.out:  simplect.c
	gcc simplect.c -lm

plot10L.o:plot.py a.out
	python plot.py 'condiciones9.txt' -0.2 0.2 -0.2 0.2 Simplectico_10

11R.out:  RK4.c
	gcc RK4.c -lm

plot11R.o:plot.py a.out
	python plot.py 'condiciones9.txt' -3 -1.5 -1 1 RK4_11
plot17R.o:ploth.py datos.dat
	python ploth.py E_RK4

11L.out:  simplect.c
	gcc simplect.c -lm

plot11L.o:plot.py a.out
	python plot.py 'condiciones9.txt' -3 -1.5 -1 1 Simplectico_11

plot17L.o:ploth.py datos.dat
	python ploth.py E_Simplect

clean:
	rm *.out *.dat *png


