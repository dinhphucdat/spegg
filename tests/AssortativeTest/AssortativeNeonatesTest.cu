#include<iostream>
#include"gtest/gtest.h"
#include"GenPhenMapTest.h"
#include<species/inds_stochastic.h>
#include<species/add_kids/neonates_class.h>
#include<species/add_kids/assortative_mating_neonates_class.h>
#include<species/add_kids/parents_class.h>

class Assortative_mating_neonate_test : public Assortative_mating_neonates {
    public:
        Assortative_mating_neonate_test(inds_stochastic* species, thrust::device_vector<int> &kids_per_mom) : 
            Assortative_mating_neonates(species, kids_per_mom)
        {}

        thrust::device_vector<int> get_mothers_chosen()
        {
            return mothers_chosen;
        }
};

class TestParents : public Parents {
    public:
        int FECUNDITY_PHENOTYPE_INDEX;
        TestParents(inds_stochastic *species) : Parents(species)
        {
            // std::cout << "before do stuff in MothParents()" << std::endl;
            FECUNDITY_PHENOTYPE_INDEX = demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
            // std::cout << "after do stuff in MothParents()" << std::endl;
            thrust::fill(kids_per_mom.begin(), kids_per_mom.end(), 2);
        }

        void determine_probability_individual_becomes_female_parent()
        {
            thrust::multiplies<float> op;
            // thrust::transform(will_reproduceF.begin(), will_reproduceF.begin() + size, phenotype[FECUNDITY_PHENOTYPE_INDEX].begin(),  probability_individual_becomes_female_parent.begin(), op);
            thrust::transform(will_reproduceF.begin(), will_reproduceF.begin() + size, phenotype[FECUNDITY_PHENOTYPE_INDEX].begin(), probability_individual_becomes_female_parent.begin(), op);
        }

        void determine_probability_individual_becomes_male_parent()
        {
            thrust::multiplies<float> op;
            thrust::transform(will_reproduceM.begin(), will_reproduceM.begin() + size, phenotype[FECUNDITY_PHENOTYPE_INDEX].begin(), probability_individual_becomes_male_parent.begin(), op);
        }

};

class TestSpecies : public inds_stochastic, public ::testing::Test {
    public:
        friend class Assortative_mating_neonate_test;
        Assortative_mating_neonate_test* neonate;
        TestParents* parents;
        TestSpecies() : TestSpecies(10, 40, 1, 4, 0) {

        }
        TestSpecies(int sizeVal, int maxSizeVal, int seedVal, int nDemes, int speciesIDVal) : 
            inds_stochastic(sizeVal, maxSizeVal, seedVal, nDemes, speciesIDVal)
        {
            // assigns every individual with an age of 0.
            thrust::fill(age.begin(), age.begin() + size, 0);
            // init demes
            initializeDemes();
            // Specify the genetics by assuming allelic values are gaussian-distributed
            for (int i = 0; i < nloci; i++)
            {
                draw_gaussian(size, 0, 0.1, fgenotype[i], gen);
                draw_gaussian(size, 0, 0.1, mgenotype[i], gen);
            }
            // std::cout << "inside moths 2" << std::endl;
            // std::cout << "at initialization size is " << size << std::endl;
            //set phenotypes
            setPhenotype(0, size);
            // std::cout << "inside moths 3. Size is " << size << std::endl;
            // To start, assign odd numbered individuals to be male, even numbered individuals to be female
            thrust::device_vector<int> twos(size);
            thrust::fill(twos.begin(), twos.end(), 2);
            thrust::transform(id.begin(), id.begin() + size, twos.begin(), sex.begin(), thrust::modulus<int>());
            demeCalculations();
            parents = new TestParents(this);
            parents->setup_parents();
            neonate = new Assortative_mating_neonate_test(this, parents->kids_per_mom);
        }

        void initializeDemes()
        {
            thrust::fill(max_deme_sizes.begin(), max_deme_sizes.begin() + Num_Demes, maxsize/Num_Demes);

            for (int i = 0; i < Num_Demes; i++)
            {
                float startSize = (float) size/Num_Demes;
                int temp1 = (int) startSize * i;
                int temp2 = (int) startSize * (i + 1);
                thrust::fill(deme.begin() + temp1, deme.begin() + temp2, i);
            }
        }
        void setPhenotype(int index, int numIndsToCalc)
        {
            for (int i = 0; i < nphen; i++)
            {
                GenotypePhenotypeMap* genPhenMap;
                genPhenMap = genPhenMap->create_genotype_phenotype_map(this, i, index, numIndsToCalc);
                genPhenMap->calculate_phenotype(this);
                delete genPhenMap;
            }
        }
        ~TestSpecies()
        {
            delete neonate;
            delete parents;
        }
};

::testing::AssertionResult expectEq(const char* expr1, const char* expr2, int expected, int actual)
{
    return ((expected == actual) ? ::testing::AssertionSuccess() : ::testing::AssertionFailure()) << "Expected = " << expected << " || Actual = " << actual;
}

TEST_F(TestSpecies, AssortativeNeonatePairsDeme) {
    neonate->inherit_genotypes_by_pair(parents->probability_individual_becomes_female_parent, 
        deme, 
        deme, 
        fgenotype,
        mgenotype);
    thrust::device_vector<int> mothers_chosen = neonate->get_mothers_chosen();
    ASSERT_NE(neonate->get_mothers_chosen().size(), 0);
    ASSERT_NE(neonate->kids_deme.size(), 0);
    ASSERT_EQ(neonate->get_mothers_chosen().size(), neonate->kids_deme.size());
    // for (int i = 0; i < 15; i++) // --> this is gonna pass the test
    for (int i = 0; i < mothers_chosen.size(); i++) {
        SCOPED_TRACE("Testing index " + std::to_string(i) + " - will be interrupted if there is an error");
        EXPECT_PRED_FORMAT2(expectEq, neonate->kids_deme[i], mothers_chosen[i]);
    }
}


class TestSpeciesNoAssort : public inds_stochastic, public ::testing::Test {
    public:
        EggsNeonates* neonate;
        TestParents* parents;
        TestSpeciesNoAssort() : TestSpeciesNoAssort(10, 40, 1, 4, 0) {

        }
        TestSpeciesNoAssort(int sizeVal, int maxSizeVal, int seedVal, int nDemes, int speciesIDVal) : 
            inds_stochastic(sizeVal, maxSizeVal, seedVal, nDemes, speciesIDVal)
        {
            // assigns every individual with an age of 0.
            thrust::fill(age.begin(), age.begin() + size, 0);
            // init demes
            initializeDemes();
            // Specify the genetics by assuming allelic values are gaussian-distributed
            for (int i = 0; i < nloci; i++)
            {
                draw_gaussian(size, 0, 0.1, fgenotype[i], gen);
                draw_gaussian(size, 0, 0.1, mgenotype[i], gen);
            }
            // std::cout << "inside moths 2" << std::endl;
            // std::cout << "at initialization size is " << size << std::endl;
            //set phenotypes
            setPhenotype(0, size);
            // std::cout << "inside moths 3. Size is " << size << std::endl;
            // To start, assign odd numbered individuals to be male, even numbered individuals to be female
            thrust::device_vector<int> twos(size);
            thrust::fill(twos.begin(), twos.end(), 2);
            thrust::transform(id.begin(), id.begin() + size, twos.begin(), sex.begin(), thrust::modulus<int>());
            demeCalculations();
            parents = new TestParents(this);
            parents->setup_parents();
            neonate = new EggsNeonates(this, parents->kids_per_mom);
        }

        void initializeDemes()
        {
            thrust::fill(max_deme_sizes.begin(), max_deme_sizes.begin() + Num_Demes, maxsize/Num_Demes);

            for (int i = 0; i < Num_Demes; i++)
            {
                float startSize = (float) size/Num_Demes;
                int temp1 = (int) startSize * i;
                int temp2 = (int) startSize * (i + 1);
                thrust::fill(deme.begin() + temp1, deme.begin() + temp2, i);
            }
        }
        void setPhenotype(int index, int numIndsToCalc)
        {
            for (int i = 0; i < nphen; i++)
            {
                GenotypePhenotypeMap* genPhenMap;
                genPhenMap = genPhenMap->create_genotype_phenotype_map(this, i, index, numIndsToCalc);
                genPhenMap->calculate_phenotype(this);
                delete genPhenMap;
            }
        }
        ~TestSpeciesNoAssort()
        {
            delete neonate;
            delete parents;
        }
};


TEST_F(TestSpeciesNoAssort, AssortativeNeonatePairsDeme) {
    neonate->inherit_genotypes(parents->probability_individual_becomes_female_parent, 
        parents->probability_individual_becomes_male_parent);
    ASSERT_NE(neonate->kids_deme.size(), 0);
    // for (int i = 0; i < 15; i++) {// --> this is gonna pass the test
    for (int i = 0; i < neonate->kids_deme.size(); i++) {
        SCOPED_TRACE("Testing index " + std::to_string(i) + " - will be interrupted if there is an error");
        EXPECT_PRED_FORMAT2(expectEq, neonate->kids_deme[i], deme[maternal_id[i + neonate->previous_pop_size]]);
    }
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
