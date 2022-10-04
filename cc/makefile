gslinc = $(shell gsl-config --cflags)
gsllib = $(shell gsl-config --libs)

libs: libk3cd.so

libk3cd.so : scatter.c
	gcc ./scatter.c -Ofast -fPIC -shared $(gslinc) $(gsllib) -o ./libk3cd.so
