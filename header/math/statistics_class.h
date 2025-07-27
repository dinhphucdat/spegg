#ifndef STATISTICS_H
#define STATISTICS_H

#include <species/inds.h>
#include <species/inds_stochastic.h>
#include <util/thrust_functors.h>
#include <util/reduce_by_key_with_zeroes.h>
#include <math/histogram.h>

#include <fstream>
#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/functional.h>

/**
 * @brief A class for obtaining summary statistics about data stored in inds
 * 
 * <b>Some features this class provides</b>
 * 1. Mean, Min, and Max phenotypes by demes
 * 2. Mean genotypes by demes
 * 3. Phenotypic and Genotypic Variance by demes
 * 4. Generate histograms based on a particular phenotype
 * 5. Calculate the quantiles based on a particular phenotype
 * 
 * Moreover, this class allows users to specify the names of the blank files to output the summary statistics and quantiles
 * 
 * @author xxx
 */
class Statistics
	{
	public:
		/**
		 * @brief Construct a new Statistics object and initialize necessary vector objects to store population statistics
		 * 
		 * @param num_demes Number of demes
		 */
		Statistics(int num_demes);
		/**
		 * @brief Construct a new Statistics object, receives output file names along with number of demes
		 * 
		 * @param num_demes number of demes
		 * @param output_file_summary_statistics name of file to output the population summary statistics
		 * @param output_file_quantiles name of file to output the quantiles, aka, histograms
		 * 
		 * @see Statistics::Statistics
		 */
		Statistics(int num_demes, const char *output_file_summary_statistics, const char *output_file_quantiles);
		/**
		 * @brief A simpler version of the constructor which should be used when only the quantiles (histograms) should be backed up
		 * 
		 * @param num_demes number of demes
		 * @param output_file_quantiles name of file to output the quantiles, aka, histograms
		 * @see Statistics::Statistics
		 */
		Statistics(int num_demes, const char *output_file_quantiles);
		/**
		 * @brief Destroy the Statistics object
		 * 
		 */
		~Statistics();
		/**
		 * @brief Calculate the mean phenotypes by demes
		 * 
		 * This function will pour the output deme-wise means into this class's private vector. Please call @link Statistics::output_results @endlink to output the data into the file(s) that has been specified in the constructor
		 * 
		 * @param individuals an @link inds @endlink object
		 * @param PHENOTYPE_TO_RECORD ID of the phenotype configured in the @c deme_config.txt
		 */
		void calculate_mean_phenotypes_by_deme(inds *individuals, int PHENOTYPE_TO_RECORD);
		/**
		 * @brief Calculate the mean genotypes by demes
		 * 
		 * This function will pour the output deme-wise means into this class's private vector. Please call @link Statistics::output_results @endlink to output the data into the file(s) that has been specified in the constructor
		 * 
		 * @param individuals an @link inds @endlink object
		 * @param GENOTYPE_TO_RECORD ID of the genotype configured in the @c deme_config.txt
		 */
		void calculate_mean_genotypes_by_deme(inds *individuals, int GENOTYPE_TO_RECORD);
		/**
		 * @brief Calculate the min and max phenotypes by demes
		 * 
		 * This function will pour the output deme-wise max and min phenotypes into this class's two private vectors. Please call @link Statistics::output_results @endlink to output the data into the file(s) that has been specified in the constructor
		 * 
		 * @param individuals an @link inds @endlink object
		 * @param PHENOTYPE_TO_RECORD ID of the phenotype configured in the @c deme_config.txt
		 */
		void calculate_min_max_phenotypes_by_deme(inds *individuals, int PHENOTYPE_TO_RECORD);
		/**
		 * @brief Calculate the phenotypic variance by demes
		 * 
		 * This function will pour the output deme-wise phenotypic variance into this class's private vector. Please call @link Statistics::output_results @endlink to output the data into the file(s) that has been specified in the constructor
		 * 
		 * @param individuals an @link inds @endlink object
		 * @param PHENOTYPE_TO_RECORD ID of the phenotype configured in the @c deme_config.txt
		 */
		void calculate_phenotypic_variance_by_deme(inds *individuals, int PHENOTYPE_TO_RECORD);
		/**
		 * @brief 
		 * 
		 * @note missing
		 * 
		 * @param individuals an @link inds @endlink object
		 */
		void calculate_genotype_variance_by_deme(inds *individuals);
		/**
		 * @brief Create a sample histograms object based on population's phenotypic information, store it in @link Statistics::histogram_by_deme @endlink
		 * 
		 * @param individuals an @link inds @endlink object
		 * @param number_of_bins how many units of information the histogram should be devided into. More specifically, this can be interpreted as the resolution of the histogram
		 * @param PHENOTYPE_TO_RECORD ID of the phenotype of interest, specified in @c deme_config.txt
		 * 
		 * @see Statistics::output_results
		 * @sa Statistics::output_histogram
		 */
		void create_sample_histograms(inds *individuals, int number_of_bins, int PHENOTYPE_TO_RECORD);
		/**
		 * @brief Calculate the quantiles based on phenotypes, store them in @link Statistics::quantiles_by_deme @endlink
		 * 
		 * @param individuals an @link inds @endlink object
		 * @param number_of_bins how many units of information the graph should be devided into. More specifically, this can be interpreted as the resolution of the graph
		 * @param PHENOTYPE_TO_RECORD ID of the phenotype of interest, specified in @c deme_config.txt
		 * 
		 * @see Statistics::output_results
		 */
		void calculate_quantiles(inds *individuals, int number_of_bins, int PHENOTYPE_TO_RECORD);
		/**
		 * @brief Output the simulation population's summary statistics into the @c output_file_summary_statistics provided at the constructor
		 * 
		 * This function should only be called after at lease one function to perform population's statistical analyses was called
		 */
		void output_results();
		/**
		 * @brief Output the simulation population's phenotypic histograms into the @c output_file_quantiles provided at the constructor
		 * 
		 * This function should only be called after at lease one function to perform population's statistical analyses, especially related to creating histograms and quantiles, was called
		 * 
		 * @param number_of_bins how many units of information the histogram should be devided into. More specifically, this can be interpreted as the resolution of the histogram
		 * 
		 * @see Statistics::create_sample_histograms
		 * @sa Statistics::calculate_quantiles
		 */
		void output_histogram(int number_of_bins); // For outputting the ECDF
		/**
		 * @brief Print mean phenotypes by demes, in sequence
		 */
		void print_mean_phenotypes_by_deme();
		/**
		 * @brief Print mean genotypes by demes, in sequence
		 * 
		 */
		void print_mean_genotypes_by_deme();

	protected:
		/// @brief file stream to output the summary statistics into
		std::ofstream summary_statistics;
		/// @brief file stream to output the histograms and quantiles into
		std::ofstream histogram_file;
		/// @deprecated
		void genotype_phenotype_map(inds **individuals);
		/// @brief total number of individuals in the population, read from the @link inds @endlink object
		int total_number_of_inds;
		/// @brief number of demes
		int number_of_demes;	
		/// @brief vector storing mean phenotypes by demes
		thrust::device_vector<float> mean_phenotypes;
		/// @brief vector storing mean genotypes by demes
		thrust::device_vector<float> mean_genotypes;
		/// @brief vector storing phenotypic variance by demes
		thrust::device_vector<float> phenotypic_variance;
		/// @brief vector storing genotypic variance by demes
		thrust::device_vector<float> genetic_variance;
		/// @brief vector storing max phenotypes by demes
		thrust::device_vector<float> max_phenotypes;
		/// @brief vector storing min phenotypes by demes
		thrust::device_vector<float> min_phenotypes;
		/// @brief vector storing size of each deme
		thrust::device_vector<float> deme_sizes;
		/// @brief dynamically allocated vector as a buffer to the output histograms
		thrust::device_vector<int> *histogram_by_deme;
		/// @brief dynamically allocated vector as a buffer to the output quantiles
		thrust::device_vector<float> *quantiles_by_deme;
	};


struct variance_elements_calculator
	{
	/* 
	Elements in the tuple.
	---------------------
	0: individual's random number value
	1: deme-wide mean
	2: deme size
	3: return value used to calculate variance
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		thrust::get<3>(t) = ((thrust::get<0>(t)-thrust::get<1>(t))*(thrust::get<0>(t)-thrust::get<1>(t)))/thrust::get<2>(t);
		}   
	};


#endif
