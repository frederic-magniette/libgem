all: libgem.so mfit.exe lfit.exe ufit.exe spfit.exe
	make -C data all
	make -C tests all

OBJS= point.o angle.o dataset.o line.o circle.o spiral.o weights.o vector.o algebra.o gaussian.o normal_gene.o  gem.o matrix.o gline.o gcircle.o object.o algos.o distrib.o gnuplot.o sdl.o graphics.o neighbouring.o

%.o : %.c libgem.h
	gcc -g -Wall -fpic -c $<

%.exe : %.c libgem.h libgem.so
	gcc -g -Wall -o $@ -L. -I. $< -lgem -lSDL -lm

libgem.so: libgem.h $(OBJS)
	gcc --shared -o libgem.so $(OBJS) 

clean:
	rm -f *~ *.o *.exe *.so *.gdat
	make -C data clean
	make -C tests clean

install:
	cp -f libgem.so /usr/local/lib
	cp -f mfit.exe lfit.exe ufit.exe /usr/local/bin

liblink:
	rm -f /usr/local/lib/libgem.so
	rm -f /usr/local/include/libgem.h
	ln -s `pwd`/libgem.so /usr/local/lib
	ln -s `pwd`/libgem.h /usr/local/include
	ldconfig
