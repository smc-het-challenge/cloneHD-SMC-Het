CC		= g++
GSL		= `gsl-config --cflags --libs`
CC_FLAGS	= -Wall -O3 -static-libgcc -static-libstdc++ -DHAVE_INLINE
LD_FLAGS	= ${GSL}
H_OBJECTS	= metrics.o 

metrics: $(H_OBJECTS)
	$(CC) $(CC_FLAGS) $(H_OBJECTS) -o run_metrics $(LD_FLAGS)
metrics.o: metrics.cpp
	$(CC) $(CC_FLAGS) -c metrics.cpp
	
clean:
	rm -f *.o
