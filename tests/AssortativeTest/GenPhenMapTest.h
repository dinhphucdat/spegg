#include <species/add_kids/genotype_phenotype_map.h>
#include <species/inds.h>

class FecundityPhenotype : public GenotypePhenotypeMap
{
    public:
        FecundityPhenotype(inds* species, int phenotypeIndex, int indexCase, int numKids) : GenotypePhenotypeMap(species, phenotypeIndex, indexCase, numKids){
            // std::cout << "inside fecundity constructor" << std::endl;
        };

        using GenotypePhenotypeMap::calculate_phenotype;
        void calculate_phenotype(inds* species);

};



struct FecundityCalculator {
    float *mapConstants;
    float *coefficient0;
    float *coefficient1;

    FecundityCalculator(float *con, float *coef0, float *coef1) : mapConstants(con), coefficient0(coef0), coefficient1(coef1){};

    /* 
    Elements in the tuple.
    ---------------------
    0: individual's deme
    1: individual's maternally inherited allelic value at their locus 0
    2: individual's paternally inherited allelic value at their locus 0
    3: individual's maternally inherited allelic value at their locus 1
    4: individual's paternally inherited allelic value at their locus 1
    5: individual's fecundity phenotypic value
    */
    template<typename tuple>
    __host__ __device__
    void operator()(tuple t)
    {
        int individualDeme = thrust::get<0>(t);

        thrust::get<5>(t) = mapConstants[individualDeme] + coefficient0[individualDeme] * (thrust::get<1>(t) + thrust::get<2>(t)) / 2 + coefficient1[individualDeme] * (thrust::get<3>(t) + thrust::get<4>(t)) / 2;

    }
};

class MortalityPhenotype : public GenotypePhenotypeMap
{
    public:
        MortalityPhenotype(inds* species, int phenotypeIndex, int indexCase, int numKids) : GenotypePhenotypeMap(species, phenotypeIndex, indexCase, numKids){};

        using GenotypePhenotypeMap::calculate_phenotype;
        void calculate_phenotype(inds* species);

};



struct MortalityCalculator {
    float* mapConstants;
    float* mapCoefficient;

    MortalityCalculator(float* con, float* coef) : mapConstants(con), mapCoefficient(coef){};

    /* 
	Elements in the tuple.
	---------------------
	0: individual's deme
	1: individual's maternally inherited allelic value at their locus 0
	2: individual's paternally inherited allelic value at their locus 0
	3: individual's mortality phenotypic value
	*/ 

    template<typename tuple>
    __host__ __device__
    void operator()(tuple t)
    {
        int individualDeme = thrust::get<0>(t);

        thrust::get<3>(t) = mapConstants[individualDeme] + mapCoefficient[individualDeme] * (thrust::get<1>(t) + thrust::get<2>(t))/2;
    }
};

class WingColorPhenotype : public GenotypePhenotypeMap
{
    public:
        WingColorPhenotype(inds* species, int phenotypeIndex, int indexCase, int numKids) : GenotypePhenotypeMap(species, phenotypeIndex, indexCase, numKids){};

        using GenotypePhenotypeMap::calculate_phenotype;
        void calculate_phenotype(inds* species);

};



struct WingColorCalculator {
    float* mapConstants;
    float* mapCoefficients;

    WingColorCalculator(float* con, float* coef) : mapConstants(con), mapCoefficients(coef){};

    /* 
	Elements in the tuple.
	---------------------
	0: individual's deme
	1: individual's maternally inherited allelic value at their locus 0
	2: individual's paternally inherited allelic value at their locus 0
	3: individual's wing color phenotypic value
	*/

    template<typename tuple>
    __host__ __device__
    void operator()(tuple t)
    {
        int individualDeme = thrust::get<0>(t);

        thrust::get<3>(t) = mapConstants[individualDeme] + mapCoefficients[individualDeme] * (thrust::get<1>(t) + thrust::get<2>(t))/2;
    }
};