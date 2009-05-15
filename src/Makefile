C_HDRS = ps_constants.h ps_errors.h ps_data.h 
C_OBJS = ps_data.o ps_data_io.o psi_shooter.o

psi_shooter: ${C_OBJS} 
	gcc -o $@ ${C_OBJS}

test_ps_data_io: ps_data.o ps_data_io.o test_ps_data_io.o
	gcc -o $@ ps_data.o ps_data_io.o test_ps_data_io.o
.c.o :
	gcc -c $<
	
clean:
	rm psi_shooter
	rm *.o