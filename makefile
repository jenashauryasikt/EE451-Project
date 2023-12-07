all: serial parallel

serial: serial_ats.cpp
	g++ serial_ats.cpp -o serial

parallel: async_read_pthreads.cpp
	g++ async_read_pthreads.cpp -o async_read_pthreads -lpthread

clean: 
	rm serial async_read_pthreads
