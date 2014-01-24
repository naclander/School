#include<stdio.h>
#include<omp.h>
#include<stdlib.h>
#include<string.h>
#include<pthread.h>

//lock used when writing to buffer
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;


int compare(const int * a, const int * b){
	if (*a > *b){
		return(1);
	}
	else if( *b > *a){
		return(-1);
	}
	else{
		return(0);
	}
}

//is_prime code required by assignment
int is_prime(int p){
	int i;
	if(p < 2){
		return(0);
	}
	i = 2;
	while(i*i <= p){
		if(p % i == 0){
			return 0;
		}
		i++;
	}
	return(1);
}

int find_primes(int max_num, int * buffer){
	int buf_index = 0;
	int i;
	#pragma omp parallel for
	for(i = 2; i <= max_num; i++){
		if(is_prime(i)){
			pthread_mutex_lock(&mutex);
			buffer[buf_index] = i;
			buf_index++;
			pthread_mutex_unlock(&mutex);
		}
	}
	return(buf_index);
}

int main(int argc, char **argv){
	if(argc < 3){
		printf("usage: hw2a MAX_NUM MAX_THREADS\n");
		return(EXIT_FAILURE);
	}
	unsigned int max_num     = atoi(argv[1]);
	unsigned int num_threads = atoi(argv[2]);
	int * buffer = malloc(max_num*sizeof(int));
	//set openmp threads
	omp_set_num_threads(num_threads);
	int num_primes = find_primes(max_num,buffer);
	int i;
	if(num_primes > 0){
		qsort(buffer,num_primes,sizeof(int),(int(*)(const void*,const void*))compare);
		printf("%d",buffer[0]);
		for(i = 1; i < num_primes; i++){
			printf(", %d", buffer[i]);
		}
		printf("\n");
	}
	return(EXIT_SUCCESS);
}
