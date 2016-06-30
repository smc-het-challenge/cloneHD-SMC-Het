CC		= g++
CC_FLAGS	= -Wall -O3 -static-libgcc -static-libstdc++ -DHAVE_INLINE
H_OBJECTS	= smchet_report.o 

metrics: $(H_OBJECTS)
	$(CC) $(CC_FLAGS) $(H_OBJECTS) -o smchet_report
metrics.o: smchet_report.cpp
	$(CC) $(CC_FLAGS) -c smchet_report.cpp
clean:
	rm -f *.o
