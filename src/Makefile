C_HDRS = ps_constants.h ps_errors.h ps_data.h
C_OBJS = ps_data.o psi_shooter.o

psi_shooter: ${C_OBJS} 
	gcc -o $@ ${C_OBJS}

.c.o :
	gcc -c $<
	
clean:
	rm psi_shooter
	rm *.o