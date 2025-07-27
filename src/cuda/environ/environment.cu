#include <environ/environment.h>

#include <curand.h>
#include <fstream>
#include <iostream>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/remove.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/scatter.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

/**
 * Intializes the space in the memory to record biotic and abiotic effects later on in the simulation and primes up the random generator.
 * @param seed_val initial seed to prime up the random generator
 * @param num_biotic_variables number of biotic factors influencing species's lives
 * @param num_abiotic_variables number of abiotic factors influencing species's lives
 * @param num_demes number of demes (individuals' clusters) in the population
 */
environment::environment(int seed_val, int num_biotic_variables, int num_abiotic_variables, int num_demes) : seed(seed_val), ndemes(num_demes), nbiotic_vars(num_biotic_variables), nabiotic_vars(num_abiotic_variables)
	{
	//Allocate biotic data vectors.
	biotic_variables = new thrust::device_vector<float>[nbiotic_vars];
	effect_of_inds_on_biotic_variable = new thrust::device_vector<float>[nbiotic_vars];

	// Resize to accord with number of demes
	for (int i = 0 ; i < nbiotic_vars ; i++) 
		{
		biotic_variables[i].resize(ndemes);
		effect_of_inds_on_biotic_variable[i].resize(ndemes);
		}

		// Initialize the first part of the feedback process 
	for (int i = 0 ; i < nbiotic_vars ; i++) 
		{
		thrust::fill(effect_of_inds_on_biotic_variable[i].begin(), effect_of_inds_on_biotic_variable[i].begin() + ndemes, 0);
		}

	// Allocate abiotic data vectors
	abiotic_variables = new thrust::device_vector<float>[nabiotic_vars];

	for (int i=0; i < nabiotic_vars; i++)
		{
		abiotic_variables[i].resize(ndemes);
		}

	//Initialize curand generator.
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
	curandSetPseudoRandomGeneratorSeed(gen, seed);
	}

/**
 * Sets up the abiotic matrix.
 * 
 * @note
 * > <ul>An abiotic matrix including three vectors, both of same size
 * >    <li> @c abiotic_variable_names list of abiotic factors' names</li>
 * >    <li> @c abiotic_variables list of lists of abiotic values, with the inner lists whose order relative to other inner lists corresponds to the abiotic names and whose size is of deme size</li>
 * >    <li> @c abiotic_variable_indices list of indices of abiotic factors' names, used for indexing the order of abiotic factors</li>
 * > </ul>
 * 
 * @param filename name of the configuration file for environment
 * 
 */
void environment::initialize_abiotic_variables(const char *filename)
        {
        Config cfg;
        try
                {
                cfg.readFile(filename);
                }
        catch(const FileIOException &fioex)
                {
                std::cerr << "No environmental config file available" << std::endl;
                }

        const Setting &root = cfg.getRoot();

        const Setting &abiotic_variable_specification = root["abiotic_variable_names"];

        int lenVal = abiotic_variable_specification.getLength();

        abiotic_variable_names.resize(lenVal);

        for (int i=0; i < lenVal; i++)
                {
                abiotic_variable_names[i] = abiotic_variable_specification[i].c_str();
                }

        // Read in the values for each abiotic variable

        const Setting &abiotic_variable_values = root["abiotic_variables"] ;

        abiotic_variables = new thrust::device_vector<float>[nabiotic_vars];

        for (int i=0; i < nabiotic_vars; i++)
                {
                abiotic_variables[i].resize(ndemes);
                }

        for (int i=0; i < nabiotic_vars; i++)
                {
                for (int j=0; j < ndemes; j++)
                        {
                        const Setting &deme_values = abiotic_variable_values[j];
                        float val = 0;
                        deme_values.lookupValue(abiotic_variable_names[i], val);
                        abiotic_variables[i][j] = val;
                        }
                }

        // specify the indices associated with each abiotic variable's name:

        for (int i=0; i < nabiotic_vars; i++)
                {
                abiotic_variable_indices[abiotic_variable_names[i]] = i;
                }
        }

/**
 * Gets the pointer to a vector of deme-wise abiotic values of one abiotic factor specified by its name
 * 
 * @param abiotic_variable_name the name of the abiotic factor whose vector of its set of values is desired
 * @return a float vector of the deme-wise values of that abiotic factor
 * @sa environment::initialize_abiotic_variables
 */
thrust::device_ptr<float> environment::get_abiotic_vector_ptr(const char *abiotic_variable_name)
	{
	return(&abiotic_variables[abiotic_variable_indices[abiotic_variable_name]][0]);
	}

/**
 * Destructor
 */
environment::~environment()
	{
	delete[] biotic_variables;
	delete[] effect_of_inds_on_biotic_variable;
	}


