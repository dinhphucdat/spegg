#ifndef GENOTYPE_PHENOTYPE_MAP
#define GENOTYPE_PHENOTYPE_MAP

#include <species/inds.h>
#include <species/add_kids/genotype_phenotype_map_parameters.h>

//virtual interface

/**
 * @brief An interface (pure virtual class) to calculate the phenotype from genotypes
 * 
 */
class GenotypePhenotypeMap
	{
	public:
/* use a factory method */
		/**
		 * @brief Create a genotype phenotype map object. This serves as a factory method to create derived classes of @link GenotypePhenotypeMap @endlink
		 * 
		 * @param species an @link inds @endlink object
		 * @param phenotype_index phenotype index
		 * @param index_case index case, which specifies the type of genotype-phenotype map to be used
		 * @param num_kids number of kids to be produced
		 * @return GenotypePhenotypeMap* 
		 */
		static GenotypePhenotypeMap *create_genotype_phenotype_map(inds *species, int phenotype_index, int index_case, int num_kids);
		/**
		 * @brief Construct a new Genotype Phenotype Map object
		 * 
		 * @param species an @link inds @endlink object
		 * @param phenotype_index phenotype index
		 * @param index_case index case, which specifies the type of genotype-phenotype map to be used
		 * @param num_kids number of kids to be produced
		 */
		GenotypePhenotypeMap(inds *species, int phenotype_index, int index_case, int num_kids)
			{
			this->phenotype_index = phenotype_index;
			this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[phenotype_index];
			this->index_case = index_case;
			this->num_kids = num_kids;
			}
		/**
		 * @brief Destroy the Genotype Phenotype Map object
		 * 
		 * @param species an @link inds @endlink object
		 */
		virtual void calculate_phenotype(inds *species)=0;

	protected:
		/// @brief pointer to the genotype-phenotype map parameters
		GenotypePhenotypeMapParameters *Parameters;
		/// @brief phenotype index
		int phenotype_index;
		/// @brief index case, which specifies the type of genotype-phenotype map to be used
		int index_case;
		/// @brief number of kids to be produced
		int num_kids;
		/// @brief random number generator
		curandGenerator_t gen;		
	};

#endif
