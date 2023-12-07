# EE451-Project

### Complie the serial code
`g++ serial_ats.cpp -o serial_ats`

### Run the serial code
`./serial_ats`

### Compile the parallel code
`g++ async_read_pthreads.cpp -o async_read_pthreads -lpthread`

### Run the parallel code
`./async_read_pthreads <NUM_OF_WORKER_BLOCKS>`
