CC				= g++
CC_FLAGS	= -Wall -O3 -static-libgcc -static-libstdc++ -DHAVE_INLINE
H_OBJECTS	= cloneHD_smchet_report.o 

metrics: $(H_OBJECTS)
	$(CC) $(CC_FLAGS) $(H_OBJECTS) -o cloneHD_smchet_report
metrics.o: cloneHD_smchet_report.cpp
	$(CC) $(CC_FLAGS) -c cloneHD_smchet_report.cpp
clean:
	rm -f *.o