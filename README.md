# ParallelMatrixMult
A simple parallel dense matrix in C/C++ using POSIX threads

This project is continously under construction.

Matrices are represented as strucuts with the number of columns and rows, an array of data and flag indicating if the matrix 
is row or column ordered.

I use POSIX threads in order to implement the parallel matrix multiplicaton. The parallelizaiton for computing A*B is 
accomplished by splitting the rows of A across the different threads specified by the user in the preprocessor 
directive NUM_THREADS. The times for the sequential and parallel implementation are reported. 

The size of matrices can be changed  by modifying the array msizes. The matrices are randomly generated by using the 
std::rand function.
