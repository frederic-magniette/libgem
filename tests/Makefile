all: isotropy.exe angres.exe perf.exe cer.exe perf_noise.exe test_algebra.exe test_gaussian.exe test_comp.exe test_gui.exe perf_algo.exe perf_algo_wc.exe 

%.exe : %.c ../libgem.h ../libgem.so
	gcc -g -Wall -fPIC -o $@ -L.. -I.. $< -lgem -lSDL -lm

clean:
	rm -f *~ *.o *.exe *.so *.gdat

install:
