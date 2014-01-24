#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<pthread.h>
#include <unistd.h>
#include<time.h>
#include<string.h>

#define FIFO_SIZE 32
//lock used when writing to buffer
pthread_mutex_t buffer_mutex = PTHREAD_MUTEX_INITIALIZER;
//lock used when incrementing num_in
pthread_mutex_t num_in_mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
	double x;
	double y;
} points;

typedef struct {
	points points_buf[FIFO_SIZE];
	int read;
	int write;
	int done;
	int num_elms;
}fifo_buf;

void init_fifo_buf(fifo_buf * buffer){
	buffer->read = 0;
	buffer->write = 0;
	buffer->done = 0;
	buffer->num_elms = 0;
}

//Add points to FIFO buffer
void add_points(points * current_points, fifo_buf * buffer){
	while(1){
		pthread_mutex_lock(&buffer_mutex);
		if(buffer->read == buffer->write && buffer->num_elms != 0){
			pthread_mutex_unlock(&buffer_mutex);
			continue;
		}
		break;
	}
	//printf("Writing points\n");
	buffer->points_buf[buffer->write].x = current_points->x;
	buffer->points_buf[buffer->write].y = current_points->y;
	buffer->write = (buffer->write + 1) % FIFO_SIZE;
	buffer->num_elms += 1;
	pthread_mutex_unlock(&buffer_mutex);
	return;
}

//Read points from FIFO buffer
void read_points(points * current_points,fifo_buf * buffer){
	while(1){
		pthread_mutex_lock(&buffer_mutex);
		if(buffer->write == buffer->read && buffer->num_elms == 0){
			pthread_mutex_unlock(&buffer_mutex);
			continue;
		}
		break;
	}	
	//printf("Reading points\n");
	current_points->x = buffer->points_buf[buffer->read].x;
	current_points->y = buffer->points_buf[buffer->read].y;
	buffer->read = (buffer->read + 1) % FIFO_SIZE;
	buffer->num_elms -= 1;
	pthread_mutex_unlock(&buffer_mutex);
	return;
}
double pi_approx(int num_points){
 	 srand (time(NULL));
	 fifo_buf shared_buf;
	 init_fifo_buf(&shared_buf);

	int num_in = 0;
	int i;
	#pragma omp parallel
	{
		#pragma omp sections
		{
			#pragma omp section //Producer thread
			{
				for(i = 0; i < num_points; i++){
					//printf("Producer for-loop\n");
					points new_points;
					new_points.x = (double)rand() / (double)RAND_MAX;
					new_points.y = (double)rand() / (double)RAND_MAX;
					add_points(&new_points, &shared_buf);
				}
				shared_buf.done = 1;
			}
			#pragma omp section //Consumer thread
			{
				while(shared_buf.done == 0 || shared_buf.num_elms > 0){
					points current_points;
					read_points(&current_points, &shared_buf);
					double r = sqrt( pow((current_points.x - 0.5),2) + 
									 pow((current_points.y - 0.5),2));
					if( r <= 0.5){
						pthread_mutex_lock(&num_in_mutex);
						num_in++;
						pthread_mutex_unlock(&num_in_mutex);
					}
				}
			}
		}//end sections
	}//end parallel
	return((double)4*num_in/num_points);
}

int main(int argc, char **argv){
	if(argc < 2){
		printf("usage: hw2c NUM_POINTS\n");
		return(EXIT_FAILURE);
	}
	int num_points = atoi(argv[1]);
	printf("%.*e\n",pi_approx(num_points));
	return(EXIT_SUCCESS);
}
