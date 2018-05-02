all : noise.exe random_line.exe random_circle.exe angle_line.exe random_spiral.exe

%.exe : %.c 
	gcc -g -Wall -fPIC -o $@ -L.. -I.. $<  -lgem -lm -lSDL

clean:
	rm -f *.exe *.gdat *~
	rm -rf  data*.ignore check_*.gp check.txt
