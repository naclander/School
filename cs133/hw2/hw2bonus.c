#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

void find_primes(int max_num, char * table){
	int i;
	int sqrt_max_num = sqrt(max_num);
	#pragma omp parallel for private(i)
	for(i = 2; i <= sqrt_max_num; i++){
		if(table[i] == 1){
			int j;
			for(j = i*i; j <= max_num; j+= i){
				table[j] = 0;
			}
		}
	}
}

void print_primes(int max_num, char * table){
	char * buffer = malloc(max_num*sizeof(char));
	sprintf(buffer,"%d",2);
	int buffer_index = 1;
	char * temp_buf = malloc(1024*sizeof(char));
	int i;
	for(i = 3; i < max_num; i++){
		if(i % 100000 == 0){
			memcpy(buffer + buffer_index,"\0",1);
			printf(buffer);
			memset(buffer,0,max_num*sizeof(char));
			buffer_index = 0;
		}
		if(table[i] == 1){
			int num_printed = sprintf(temp_buf,", %d",i);
			memcpy(buffer + buffer_index, temp_buf, num_printed);
			buffer_index += num_printed;
		}
	}

	printf(buffer);
	printf("\n");
}

int main(int argc, char **argv){
	if(argc < 3){
		printf("usage: MAX_NUM NUM_THREADS\n");
		return(EXIT_FAILURE);
	}
	int max_num = atoi(argv[1]);
	int num_threads = atoi(argv[2]);
	omp_set_num_threads(num_threads);
	char * table = malloc(max_num * sizeof(char));
	memset(table,1,max_num*sizeof(char));
	int i;
	find_primes(max_num, table);
	print_primes(max_num, table);
	return(EXIT_SUCCESS);
}
