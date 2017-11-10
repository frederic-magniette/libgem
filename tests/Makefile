all: isotropy.exe angular_resolution.exe perf.exe cer.exe perf_noise.exe test_algebra.exe test_gaussian.exe test_comp.exe test_gui.exe

%.exe : %.c ../libgem.h ../libgem.so
	gcc -g -Wall -fPIC -o $@ -L.. -I.. -lgem -lSDL $<

clean:
	rm -f *~ *.o *.exe *.so *.gdat

install:
