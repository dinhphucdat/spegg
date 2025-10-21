#ifndef FOOTIMER2_H
#define FOOTIMER2_H

#include <sys/time.h>
#include <cuda.h>
#include <stdio.h>

/** @{ @name GPU Performance Diagnosis */

/**
 * @brief A function to diagnose the memory usage of the GPU
 * 
 */
void memory_diagnostic();

/**
 * @brief A class to measure elapsed time for code execution
 * 
 */
struct footimer2
	{
	public:
		/**
		 * @brief Start the timer
		 * 
		 */
		void start();
		/**
		 * @brief Stop the timer
		 * 
		 */
		void stop();
		/**
		 * @brief Print the elapsed time in seconds
		 * 
		 */
		void printTime();
		/**
		 * @brief Print the elapsed time in microseconds
		 * 
		 */
		void uprintTime();
		/**
		 * @brief Get the elapsed time in seconds
		 * 
		 * @return int 
		 */
		int getElapsed();
	private:
		struct timeval start_time;
		struct timeval stop_time;
	};

/** @} */

#endif
