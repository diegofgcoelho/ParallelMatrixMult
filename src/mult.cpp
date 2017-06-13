/*
 * mult.cpp
 *
 *  Created on: May 30, 2017
 *      Author: Diego Coelho, PhD Candidate, UofC.
 */


#include <iostream>
#include <ctime>
#include <string>
#include <cstring>
#include <pthread.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>

#define SUCCESS 0
#define FAIL -1
#define FILE_FAIL -2
#define REP 1000 //The number of replicates
#define ROW_ORDERED true
#define COL_ORDERED false

#define NUM_THREADS 4

using namespace std;

//Defining my own matrix type
typedef struct {
	//Array of data and must have size nrows*ncols
	double *data;
	//The number of rows and cols
	unsigned int nrows, ncols;
	//Flag indicating format: true if it is row ordered, and false if it is col ordered
	bool format;
} Matrixd;

//Defining struct type containing data to be sent to threads
typedef struct {
	Matrixd *A, *B, *C;
	unsigned int startRow, endRow;
	pthread_mutex_t* mutex;
} MatInfo;

int mult(Matrixd* A, Matrixd* B, Matrixd* C);
int mult_thread(Matrixd* A, Matrixd* B, Matrixd* C, unsigned int nt);
void *mult_thread_p(void* data);
double errnorm(Matrixd* A, Matrixd* B);
double randomd();
timespec diff_time(timespec init, timespec end);
double get_millisecs(timespec time);

int main(int argc, char* argv[]){

	//Time variables
	timespec time_init, time_end;
	double time_usual=0.0, time_thread=0.0;

	unsigned int msizes[] = {50, 100, 200, 300};
	unsigned int nmsizes = 4;

	//Running over different sizes
	for(unsigned int i = 0; i < nmsizes; i++){
		cout << "Entering size " << msizes[i] << endl;
		//Running over all the replicates for the size msizes[i]
		for(unsigned int j = 0; j < REP; j++){
			//Defining matrices to be used for simulation
			Matrixd A, B, Cseries, Cpara;
			A.nrows = msizes[i];
			A.ncols = msizes[i];
			A.format = ROW_ORDERED;
			A.data = new double[msizes[i]*msizes[i]];
			generate_n(A.data, msizes[i]*msizes[i], randomd);
			//Defining matrices to be used for simulation
			B.nrows = msizes[i];
			B.ncols = msizes[i];
			B.format = ROW_ORDERED;
			B.data = new double[msizes[i]*msizes[i]];
			generate_n(B.data, msizes[i]*msizes[i], randomd);
			//Defining matrices to be used for simulation
			Cseries.nrows = msizes[i];
			Cseries.ncols = msizes[i];
			Cseries.format = ROW_ORDERED;
			Cseries.data = new double[msizes[i]*msizes[i]];
			fill_n(Cseries.data, msizes[i]*msizes[i], 0.0);
			Cpara.nrows = msizes[i];
			Cpara.ncols = msizes[i];
			Cpara.format = ROW_ORDERED;
			Cpara.data = new double[msizes[i]*msizes[i]];
			fill_n(Cpara.data, msizes[i]*msizes[i], 0.0);

			//Calling series multiplication
			clock_gettime(CLOCK_MONOTONIC, &time_init);
			mult(&A, &B, &Cseries);
			clock_gettime(CLOCK_MONOTONIC, &time_end);
			time_usual += get_millisecs(diff_time(time_init, time_end));
			//Calling multiplication inside POSIX thread
			clock_gettime(CLOCK_MONOTONIC, &time_init);
			mult_thread(&A, &B, &Cpara, NUM_THREADS);
			clock_gettime(CLOCK_MONOTONIC, &time_end);
			time_thread += get_millisecs(diff_time(time_init, time_end));

			if(errnorm(&Cseries, &Cpara) > 1e-4){
				cout << "Possible mismatch." << endl;
			}

			//Deleting the elements that were created
			delete [] A.data;
			delete [] B.data;
			delete [] Cseries.data;
			delete [] Cpara.data;

		}
	}

	cout << "The times are: " << endl;
	cout << "Usual: " << time_usual/REP << "." << endl;
	cout << "Thread: " << time_thread/REP << "." << endl;
	return SUCCESS;
}


int mult(Matrixd* A, Matrixd* B, Matrixd* C) {
	/*
	 * Input:
	 * A and B are pointers for a matrix struct
	 * Output:
	 *
	 * C is a pointer for a matrix struct
	 * Description:
	 * this function fills the matrix pointed by C with the
	 * result f the matrix product A*B (in this order). The memory for
	 * C must be allocated before the calling of this function and must have
	 * all the elements equals to 0.0.
	 */

	if(!(A->format & B->format & C->format)){
		cout << "Error: all the matrices must be row ordered. Column ordered is still not supported." << endl;
		return FAIL;
	}

	for(unsigned int i = 0; i < A->nrows; i++){

		//Getting the ith row of A and C
		double *arow = &(A->data[i*A->nrows]);
		double *brow = &(B->data[i*B->nrows]);

		for(unsigned int j = 0; j < A->ncols; j++){

			double *crow = &(C->data[j*C->nrows]);

			for(unsigned int k = 0; k < B->ncols; k++){
				crow[k] += arow[j]*brow[k];
			}
		}
	}
	return SUCCESS;
}

int mult_thread(Matrixd* A, Matrixd* B, Matrixd* C, unsigned int nt) {

	//Defining mutex to control the update of the output matrix
	pthread_mutex_t mutex;

	MatInfo* datarray = new MatInfo[nt];
	pthread_t *mythreads = new pthread_t[nt];
	pthread_attr_t myattr;
	pthread_attr_init(&myattr);
	pthread_attr_setdetachstate(&myattr, PTHREAD_CREATE_JOINABLE);

	unsigned int strideRows = (unsigned int) (A->nrows/nt);
	//Note the for interval: the last iteration have to be done in separate
	for(unsigned int i = 0; i < nt-1; i++){
		MatInfo* mydata = &datarray[i];
		mydata->A = A;
		mydata->B = B;
		mydata->C = C;
		mydata->startRow = i*strideRows;
		mydata->endRow = (i+1)*strideRows-1;
		mydata->mutex = &mutex;

		int rcreate = pthread_create(&mythreads[i], &myattr, mult_thread_p, (void*) mydata);
		if (rcreate){
			cout << "ERROR: return code from pthread_create() is " << rcreate << endl;
			return FAIL;
		}
	}
	//Last iteration done in separate
	MatInfo* mydata = &datarray[nt-1];
	mydata->A = A;
	mydata->B = B;
	mydata->C = C;
	mydata->startRow = (nt-1)*strideRows;
	mydata->endRow = A->nrows-1;
	mydata->mutex = &mutex;
	int rcreate = pthread_create(&mythreads[nt-1], &myattr, mult_thread_p, (void*) mydata);
	if (rcreate){
		cout << "ERROR: return code from pthread_create() is " << rcreate << endl;
		return FAIL;
	}


	for(unsigned int i = 0; i < nt; i++){
		void *res = NULL;
		rcreate = pthread_join(mythreads[i], &res);
		if(rcreate != SUCCESS){
			cout << "ERROR: thread "<< i << " returned with code " << rcreate << endl;
			return FAIL;
		}
	}

	//Freeing the thread attribute and mutex
	pthread_attr_destroy(&myattr);
	pthread_mutex_destroy(&mutex);
	//Deleting the created arrays
	delete [] mythreads;
	delete [] datarray;

	return SUCCESS;
}

//This version of mult_thread_p is working, but I am trying to optimize
//it in the following piece of code after this commented block
//void *mult_thread_p(void* data){
//	MatInfo* mydata = (MatInfo*) data;
//
//	Matrixd *A = mydata->A;
//	Matrixd *B = mydata->B;
//	Matrixd *C = mydata->C;
//	unsigned int startRow = mydata->startRow;
//	unsigned int endRow = mydata->endRow;
//	//pthread_mutex_t *mutex = mydata->mutex;//not used
//
//
//	//It will store a copy of the rows of A
//	unsigned int sizeTempA = (endRow-startRow+1)*A->ncols;
//	double *tempA = new double[sizeTempA];
//	//It will store the data of B
//	unsigned int sizeTempB = B->nrows*B->ncols;
//	double *tempB = new double[sizeTempB];
//	//It will be copied to the rows of C
//	unsigned int sizeTempC = (endRow-startRow+1)*C->ncols;
//	double *tempC = new double[sizeTempC];
//
//	//Copying the appropriate elements
//	memcpy(tempA, &A->data[startRow*A->ncols], (endRow-startRow+1)*A->ncols);
//	memcpy(tempB, &B->data[0], B->nrows*B->ncols);
//	memcpy(tempC, &C->data[startRow*C->ncols], (endRow-startRow+1)*C->ncols);
//
//	for(unsigned int i = 0; i < (endRow-startRow+1); i++){
//		//Getting the ith row of A and C
//		double *arow = &tempA[i*A->ncols];
//		double *crow = &tempC[i*C->ncols];
//
//		for(unsigned int j = 0; j < B->ncols; j++){
//			double *brow = &tempB[j*B->nrows];
//
//
//			for(unsigned int k = 0; k < B->ncols; k++){
//				crow[k] += arow[j]*brow[k];
//			}
//		}
//	}
//
//	//Deleting the elements created here
//	delete [] tempA;
//	delete [] tempB;
//	delete [] tempC;
//
//	int* res = new int;
//	*res = SUCCESS;
//	pthread_exit((void*) res);
//}

void *mult_thread_p(void* data){
	MatInfo* mydata = (MatInfo*) data;

	Matrixd *A = mydata->A;
	Matrixd *B = mydata->B;
	Matrixd *C = mydata->C;
	unsigned int startRow = mydata->startRow;
	unsigned int endRow = mydata->endRow;
	//pthread_mutex_t *mutex = mydata->mutex;//not used


	//It will store a copy of the rows of A
	unsigned int sizeTempA = (endRow-startRow+1)*A->ncols;
	double *tempA = new double[sizeTempA];
	//It will store the data of B
	unsigned int sizeTempB = B->nrows*B->ncols;
	double *tempB = new double[sizeTempB];
	//It will be copied to the rows of C
	unsigned int sizeTempC = (endRow-startRow+1)*C->ncols;
	double *tempC = new double[sizeTempC];

	//Copying the appropriate elements
	tempA = &A->data[startRow*A->ncols];
	tempB = &B->data[0];
	tempC = &C->data[startRow*C->ncols];

	for(unsigned int i = 0; i < (endRow-startRow+1); i++){
		//Getting the ith row of A and C
		double *arow = &tempA[i*A->ncols];
		double *crow = &tempC[i*C->ncols];

		for(unsigned int j = 0; j < B->ncols; j++){
			double *brow = &tempB[j*B->nrows];


			for(unsigned int k = 0; k < B->ncols; k++){
				crow[k] += arow[j]*brow[k];
			}
		}
	}

	int* res = new int;
	*res = SUCCESS;
	pthread_exit((void*) res);
}


double errnorm(Matrixd* A, Matrixd* B){
	/*
	 * Input:
	 * A and B are pointers for Matrixd struct
	 * Output:
	 * the output is the double fronorm and it represents the Frobenius norm of the difference of A-B.
	 * Description:
	 * It computes the Frobenius nor of the difference A-B
	 */

	//Sanity Check
	if((A->nrows != B->nrows) || (A->ncols != B->nrows)){
		cout << "Error: the input matrices must have the same dimension." << endl;
		return FAIL;
	}

	double fronorm = 0.0;

	for(unsigned int i = 0; i < A->nrows*A->ncols; i++) fronorm += pow(A->data[i]-B->data[i], 2);

	return sqrt(fronorm);
}

double randomd(){
	/*
	 * Input:
	 * Output:
	 * double in the interval 01 from the uniform distribution
	 * Description:
	 * It returns a double in the interval 0, 1 from the uniform distribution using the std::rand function.
	 * This is just  wrapper for the generate_n function to fill the matrices.
	 */
	return (double)(rand()/RAND_MAX);
}

timespec diff_time(timespec init, timespec end){
	/*
	 * Input:
	 * init and end are timespec structures that represents the starting and ending time of the measure
	 * that we are insterested in computing the difference
	 * Output:
	 * out_time is a timespec structure representing the time elapsed from the sart to end
	 * Description:
	 * this function computes the difference between start and end time, which means the time elapsed between the beginning
	 * of the measurement and the end
	 */

	//Output variable
	timespec out_time;

	if((end.tv_nsec-init.tv_nsec) < 0){
		out_time.tv_sec = end.tv_sec-init.tv_sec-1;
		out_time.tv_nsec = 1e9 + end.tv_nsec - init.tv_nsec;
	} else {
		out_time.tv_sec = end.tv_sec - init.tv_sec;
		out_time.tv_nsec = end.tv_nsec - init.tv_nsec;
	}

	return out_time;
}

double get_millisecs(timespec time){
	/*
	 * Input:
	 * time is a timespec structure
	 * Output:
	 * out_time is a double representins the time in milliseconds
	 * Description:
	 * this function converts the time in the timespec structure time to milliseconds measures
	 */

	double out_time = static_cast<double>(time.tv_sec)*1000+static_cast<double>(time.tv_nsec)/1000000;
	return out_time;
}
