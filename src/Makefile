C_HDRS = ps_constants.h ps_errors.h ps_data.h ps_data_io.h ps_list.h
C_OBJS = ps_data.o ps_data_io.o psi_shooter.o ps_list.o
CFLAGS = -Wall
psi_shooter: ${C_OBJS} 
	gcc ${CFLAGS} -o $@ ${C_OBJS}

test_ps_data_io: ps_data.o ps_data_io.o test_ps_data_io.o
	gcc ${CFLAGS} -o $@ ps_data.o ps_data_io.o test_ps_data_io.o

test_ps_list: ps_data.o ps_list.o ps_data_io.o test_ps_list.o
	gcc ${CFLAGS} -o $@ ps_data.o ps_list.o ps_data_io.o test_ps_list.o

.c.o :
	gcc ${CFLAGS} -c $<
	
clean:
	rm psi_shooter
	rm *.o