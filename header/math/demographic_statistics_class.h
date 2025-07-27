#ifndef DEMOGRAPHIC_STATISTICS_H
#define DEMOGRAPHIC_STATISTICS_H

#include <math/statistics_class.h>
#include <species/inds.h>

#include <thrust/functional.h>
/**
 * @brief Functions the same way as its @link Statistics @endlink parent class, except it also analyzes demographic data
 * @author xxx
 * @see Statistics
 */
class DemographicStatistics : public Statistics
	{
	public:
		/**
		 * @brief Construct a new @c DemographicStatistics object, receiving names of the files to output summary statistics and histograms
		 * 
		 * @param num_demes number of demes
		 * @param output_file_summary_statistics file to output summary statistics
		 * @param output_file_histograms file to output histograms
		 */
		DemographicStatistics(int num_demes, const char *output_file_summary_statistics, const char *output_file_histograms);
		/**
		 * @brief A simpler version of the constructor, this time it only receives the file to output summary statistics into
		 * 
		 * @param num_demes number of demes
		 * @param output_file_summary_statistics file to output summary statistics
		 */
		DemographicStatistics(int num_demes, const char *output_file_summary_statistics);
		/**
		 * @brief Destroy the Demographic Statistics object and its associated dynamically allocated vectors
		 * 
		 */
		~DemographicStatistics();

		/**
		 * @brief Calculate the deme sizes
		 * 
		 * This function takes the deme sizes from the component vector inside the @link inds @endlink object which stores the deme sizes and then pours the copied deme sizes into its buffer vector
		 * 
		 * @param species an @link inds @endlink object
		 */
		void calculate_deme_sizes(inds *species);
		/**
		 * @brief output the deme sizes and sex ratios in the summary statistics file
		 * @see DemographicStatistics::record_deme_sizes
		 * @sa DemographicStatistics::calculate_sex_ratios
		 */
		void record_demographic_data();
		/**
		 * @brief output the deme abundance (by deme) into the @c output_file_summary_statistics file specified in the constructor
		 * 
		 * @note call the @link DemographicStatistics::calculate_deme_sizes @endlink first
		 * 
		 */
		void record_deme_sizes();
		/**
		 * @brief Calculate the age distribution by demes
		 * 
		 * This function outputs the age distribution histogram into this class's buffer vector, @link Statistics::histogram_by_deme @endlink
		 * 
		 * 
		 * @param species an @link inds @endlink object
		 * @param number_of_bins histogram's resolution
		 * 
		 * @see Statistics::output_results
		 * @sa Statistics::output_histogram
		 */
		void calculate_age_distribution(inds *species, int number_of_bins);
		/**
		 * @brief Calculate the sex ratios by deme
		 * 
		 * This function also outputs the deme-wise sex ratios into this class's buffer vector
		 * 
		 * @param species an @link inds @endlink object
		 * 
		 * @see Statistics::output_results
		 */
		void calculate_sex_ratios(inds *species);

	protected: 
		/// @brief number of demes
		int number_of_demes;
		/// @brief vector as a buffer to the deme sizes
		thrust::device_vector<int> deme_abundances;
		/// @brief vector as a buffer to the deme-wise sex ratios
		thrust::device_vector<float> sex_ratios;
		/// @deprecated
		void print_header();
	
	};

#endif
