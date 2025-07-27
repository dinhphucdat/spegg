#include<iostream>
#include"gtest/gtest.h"
#include"gmock/gmock.h"
#include"GenPhenMapTest.h"
#include<species/inds_stochastic.h>
#include<species/add_kids/neonates_class.h>
#include<species/add_kids/assortative_mating_neonates_class.h>
#include<species/add_kids/parents_class.h>
#include<species/add_kids/assortative_mating_parents_class.h>

class TestAssortParents : public Assortative_mating_parents {
    public:
        int FECUNDITY_PHENOTYPE_INDEX;
        TestAssortParents(inds_stochastic *species) : Assortative_mating_parents(species)
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

class Assortative_mating_neonate_test : public Assortative_mating_neonates {
    public:
        Assortative_mating_neonate_test(inds_stochastic* species, thrust::device_vector<int> kids_per_mom) : 
            Assortative_mating_neonates(species, kids_per_mom)
        {}

        thrust::device_vector<int> get_mothers_chosen()
        {
            return mothers_chosen;
        }
        thrust::device_vector<int> get_fathers_chosen()
        {
            return fathers_chosen;
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

void distributeIndvsIntoDemes(int demeSize[], int total, int numDemes)
{
    int remainder = total % numDemes;
    for (int i = 0; i < numDemes; i++)
    {
        demeSize[i] = (int) (total / numDemes) + (remainder > 0 ? 1 : 0);
        remainder--;
    }
}

struct isOdd
{
    __host__ __device__
    bool operator()(const int& x)
    {
        return x % 2 != 0;
    }
};

class TestSpecies : public inds_stochastic, public ::testing::Test {
    public:
        friend class Assortative_mating_neonate_test;
        Assortative_mating_neonate_test* neonate;
        TestAssortParents* parents;
        TestSpecies() : TestSpecies(16, 50, 1, 4, 0) {

        }
        TestSpecies(int sizeVal, int maxSizeVal, int seedVal, int nDemes, int speciesIDVal) : 
            inds_stochastic(sizeVal, maxSizeVal, seedVal, nDemes, speciesIDVal)
        {
            // assigns every individual with an age of 0.
            thrust::fill(age.begin(), age.begin() + size, 0);
            // init demes
            initializeDemes();
            demeCalculations();
            // Specify the genetics by assuming allelic values are gaussian-distributed
            for (int i = 0; i < nloci; i++)
            {
                draw_gaussian(size, 0, 0.1, fgenotype[i], gen);
                draw_gaussian(size, 0, 0.1, mgenotype[i], gen);
            }
            // std::cout << "inside moths 2" << std::endl;
            // std::cout << "at initialization size is " << size << std::endl;
            // std::cout << "inside moths 3. Size is " << size << std::endl;
            // To start, assign odd numbered individuals to be male, even numbered individuals to be female
            thrust::device_vector<int> twos(size);
            thrust::fill(twos.begin(), twos.end(), 2);
            thrust::counting_iterator<int> count(0);
            thrust::transform(count, count + size, twos.begin(), sex.begin(), thrust::modulus<int>());
            //set phenotypes
            setPhenotype(0, size);
            parents = new TestAssortParents(this);
            parents->setup_parents();
            parents->determine_parent_pair_probability();
            neonate = new Assortative_mating_neonate_test(this, parents->kids_per_mom);
        }

        void SetUp() override 
        {
            // First assertion: make sure that no demes have 0 members
            for (int i = 0; i < demeParameters->check_number_of_demes(); i++)
            {
                ASSERT_NE(deme_sizes[i], 0) << "Index " << i << " is 0";
            }
            for (int i = 0; i < demeParameters->check_number_of_demes(); i++)
            {
                ASSERT_EQ(deme_sizes[i], 4) << "Index " << i << " is not 4";
            }
            // Second assertion: make sure the system has correct number of pair evaluations
            thrust::device_vector<int> product_deme(deme_sizes);
            thrust::transform(parents->reproductive_females_per_deme.begin(), 
                            parents->reproductive_females_per_deme.end(), 
                            parents->reproductive_males_per_deme.begin(), 
                            product_deme.begin(), 
                            thrust::multiplies<int>());
            int num_pair_evals = thrust::reduce(product_deme.begin(), product_deme.end());
            ASSERT_EQ(parents->pair_demes.size(), num_pair_evals) << "Wrong number of pair evaluations || " << "Expected = " << num_pair_evals << " || Actual = " << parents->pair_demes.size();
            // Third assertion, make sure that the female_lists and male_lists were instantiated correctly
            thrust::host_vector<int> expected_females_list = {0,0,2,2,4,4,6,6,8,8,10,10,12,12,14,14};
            thrust::host_vector<int> expected_males_list   = {1,3,1,3,5,7,5,7,9,11,9,11,13,15,13,15};
            thrust::host_vector<int> females_list = parents->females_list;
            thrust::host_vector<int> males_list = parents->males_list;
            EXPECT_THAT(females_list, ::testing::ContainerEq(expected_females_list)) << ", in terms of females_list's instantiation.";
            EXPECT_THAT(males_list, ::testing::ContainerEq(expected_males_list)) << ", in terms of males_list's instantiation.";
            // Fourth assertion, make sure that every pair evaluation gets similar chance based on common assort phenotype
            thrust::device_vector<float> expected_pair_probs_dev(num_pair_evals);
            thrust::gather(parents->females_list.begin(), 
                            parents->females_list.end(), 
                            parents->probability_individual_becomes_female_parent.begin(), 
                            expected_pair_probs_dev.begin());
            thrust::device_vector<float> oh_pt_two(num_pair_evals);
            thrust::fill(oh_pt_two.begin(), oh_pt_two.end(), 0.2f);
            thrust::device_vector<float> scaled_pair_prob(num_pair_evals);
            thrust::transform(expected_pair_probs_dev.begin(), 
                                expected_pair_probs_dev.end(), 
                                oh_pt_two.begin(), 
                                scaled_pair_prob.begin(), 
                                thrust::multiplies<float>());
            thrust::transform(expected_pair_probs_dev.begin(), 
                expected_pair_probs_dev.end(), 
                scaled_pair_prob.begin(), 
                expected_pair_probs_dev.begin(), 
                thrust::plus<float>());
            thrust::host_vector<float> expected_pair_probs = expected_pair_probs_dev;
            thrust::host_vector<float> pair_probs = parents->parental_pair_probability;
            ASSERT_NE(pair_probs.size(), 0) << "Probabilities array must not be empty!";
            EXPECT_THAT(pair_probs, ::testing::ContainerEq(expected_pair_probs)) << ", in terms of pair probs when traits are fixed among individuals.";
        }

        void initializeDemes()
        {
            int demeSizeArr[Num_Demes];
            distributeIndvsIntoDemes(demeSizeArr, size, Num_Demes);
            thrust::fill(max_deme_sizes.begin(), max_deme_sizes.begin() + Num_Demes, maxsize/Num_Demes);
            int sum_so_far = 0;
            for (int i = 0; i < Num_Demes; i++)
            {
                int start = sum_so_far;
                int end = start + demeSizeArr[i];
                thrust::fill(deme.begin() + start, deme.begin() + end, i);
                sum_so_far += demeSizeArr[i];
            }
        }
        void setPhenotype(int index, int numIndsToCalc)
        {
            for (int i = 0; i < nphen; i++)
            {
                GenotypePhenotypeMap* genPhenMap;
                genPhenMap = genPhenMap->create_genotype_phenotype_map(this, i, index, numIndsToCalc);
                genPhenMap->calculate_phenotype(this);
                // this line is to test when everybody has the same phenotype, the prob should be one
                thrust::fill(phenotype[2].begin(), phenotype[2].begin() + size, 1.0f);
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
    neonate->inherit_genotypes_by_pair(
        parents->parental_pair_probability, 
        parents->males_list, 
        parents->females_list, 
        parents->pair_demes
    );
    // test if kids are from same demes as their moms
    thrust::device_vector<int> mothers_chosen = neonate->get_mothers_chosen();
    ASSERT_NE(neonate->get_mothers_chosen().size(), 0) << "no mothers chosen, which is an unexpected behavior.";
    ASSERT_NE(neonate->kids_deme.size(), 0) << "kids_deme was so poorly constructed and handled.";
    ASSERT_EQ(neonate->get_mothers_chosen().size(), neonate->kids_deme.size()) << "Number of chosen moms must be equal to number of kids.";
    thrust::device_vector<int> mothers_chosen_demes_dev(mothers_chosen.size());
    thrust::gather(mothers_chosen.begin(), mothers_chosen.end(), deme.begin(), mothers_chosen_demes_dev.begin());
    thrust::host_vector<int> kids_deme = neonate->kids_deme;
    thrust::host_vector<int> mothers_chosen_demes = mothers_chosen_demes_dev;
    EXPECT_THAT(mothers_chosen_demes, ::testing::ContainerEq(kids_deme)) << ", in terms of the required equality between mothers' demes and kids' demes.";
    // test if kids are from same demes as their dads
    thrust::device_vector<int> fathers_chosen = neonate->get_fathers_chosen();
    ASSERT_NE(neonate->get_fathers_chosen().size(), 0) << "no fathers chosen, which is an unexpected behavior.";
    ASSERT_EQ(neonate->get_fathers_chosen().size(), neonate->kids_deme.size()) << "Number of chosen dads must be equal to number of kids.";
    thrust::device_vector<int> fathers_chosen_demes_dev(fathers_chosen.size());
    thrust::gather(fathers_chosen.begin(), fathers_chosen.end(), deme.begin(), fathers_chosen_demes_dev.begin());
    thrust::host_vector<int> fathers_chosen_demes = fathers_chosen_demes_dev;
    EXPECT_THAT(fathers_chosen_demes, ::testing::ContainerEq(kids_deme)) << ", in terms of the required equality between fathers' demes and kids' demes.";
}

class TestSpeciesAssortNoChance : public inds_stochastic, public ::testing::Test {
    public:
        friend class Assortative_mating_neonate_test;
        Assortative_mating_neonate_test* neonate;
        TestAssortParents* parents;
        TestSpeciesAssortNoChance() : TestSpeciesAssortNoChance(16, 50, 1, 4, 0) {

        }
        TestSpeciesAssortNoChance(int sizeVal, int maxSizeVal, int seedVal, int nDemes, int speciesIDVal) : 
            inds_stochastic(sizeVal, maxSizeVal, seedVal, nDemes, speciesIDVal)
        {
            // assigns every individual with an age of 0.
            thrust::fill(age.begin(), age.begin() + size, 0);
            // init demes
            initializeDemes();
            demeCalculations();
            // Specify the genetics by assuming allelic values are gaussian-distributed
            for (int i = 0; i < nloci; i++)
            {
                draw_gaussian(size, 0, 0.1, fgenotype[i], gen);
                draw_gaussian(size, 0, 0.1, mgenotype[i], gen);
            }
            // std::cout << "inside moths 2" << std::endl;
            // std::cout << "at initialization size is " << size << std::endl;
            // std::cout << "inside moths 3. Size is " << size << std::endl;
            // To start, assign odd numbered individuals to be male, even numbered individuals to be female
            thrust::device_vector<int> twos(size);
            thrust::fill(twos.begin(), twos.end(), 2);
            thrust::counting_iterator<int> count(0);
            thrust::transform(count, count + size, twos.begin(), sex.begin(), thrust::modulus<int>());
            //set phenotypes
            setPhenotype(0, size);
            parents = new TestAssortParents(this);
            parents->setup_parents();
            parents->determine_parent_pair_probability();
            neonate = new Assortative_mating_neonate_test(this, parents->kids_per_mom);
        }

        void SetUp() override 
        {
            // First assertion: make sure that no demes have 0 members
            for (int i = 0; i < demeParameters->check_number_of_demes(); i++)
            {
                ASSERT_NE(deme_sizes[i], 0) << "Index " << i << " is 0";
            }
            for (int i = 0; i < demeParameters->check_number_of_demes(); i++)
            {
                ASSERT_EQ(deme_sizes[i], 4) << "Index " << i << " is not 4";
            }
            // Second assertion: make sure the system has correct number of pair evaluations
            thrust::device_vector<int> product_deme(deme_sizes);
            thrust::transform(parents->reproductive_females_per_deme.begin(), 
                            parents->reproductive_females_per_deme.end(), 
                            parents->reproductive_males_per_deme.begin(), 
                            product_deme.begin(), 
                            thrust::multiplies<int>());
            int num_pair_evals = thrust::reduce(product_deme.begin(), product_deme.end());
            ASSERT_EQ(parents->pair_demes.size(), num_pair_evals) << "Wrong number of pair evaluations || " << "Expected = " << num_pair_evals << " || Actual = " << parents->pair_demes.size();
            // Third assertion, make sure that the female_lists and male_lists were instantiated correctly
            thrust::host_vector<int> expected_females_list = {0,0,2,2,4,4,6,6,8,8,10,10,12,12,14,14};
            thrust::host_vector<int> expected_males_list   = {1,3,1,3,5,7,5,7,9,11,9,11,13,15,13,15};
            thrust::host_vector<int> females_list = parents->females_list;
            thrust::host_vector<int> males_list = parents->males_list;
            EXPECT_THAT(females_list, ::testing::ContainerEq(expected_females_list)) << ", in terms of females_list's instantiation.";
            EXPECT_THAT(males_list, ::testing::ContainerEq(expected_males_list)) << ", in terms of males_list's instantiation.";
            // Fourth assertion, make sure that every pair evaluation gets no chance based on extreme assort phenotypes
            thrust::device_vector<float> expected_pair_probs_dev(num_pair_evals);
            thrust::gather(parents->females_list.begin(), 
                            parents->females_list.end(), 
                            parents->probability_individual_becomes_female_parent.begin(), 
                            expected_pair_probs_dev.begin());
            thrust::device_vector<float> oh_pt_two(num_pair_evals);
            thrust::fill(oh_pt_two.begin(), oh_pt_two.end(), 0.2f);
            thrust::transform(expected_pair_probs_dev.begin(), 
                                expected_pair_probs_dev.end(), 
                                oh_pt_two.begin(), 
                                expected_pair_probs_dev.begin(), 
                                thrust::multiplies<float>());
            thrust::host_vector<float> expected_pair_probs = expected_pair_probs_dev;
            thrust::host_vector<float> pair_probs = parents->parental_pair_probability;
            ASSERT_NE(pair_probs.size(), 0) << "Probabilities array must not be empty!";
            EXPECT_THAT(pair_probs, ::testing::ContainerEq(expected_pair_probs)) << ", in terms of pair probs when traits are extreme among individuals.";
        }

        void initializeDemes()
        {
            int demeSizeArr[Num_Demes];
            distributeIndvsIntoDemes(demeSizeArr, size, Num_Demes);
            thrust::fill(max_deme_sizes.begin(), max_deme_sizes.begin() + Num_Demes, maxsize/Num_Demes);
            int sum_so_far = 0;
            for (int i = 0; i < Num_Demes; i++)
            {
                int start = sum_so_far;
                int end = start + demeSizeArr[i];
                thrust::fill(deme.begin() + start, deme.begin() + end, i);
                sum_so_far += demeSizeArr[i];
            }
        }
        void setPhenotype(int index, int numIndsToCalc)
        {
            for (int i = 0; i < nphen; i++)
            {
                GenotypePhenotypeMap* genPhenMap;
                genPhenMap = genPhenMap->create_genotype_phenotype_map(this, i, index, numIndsToCalc);
                genPhenMap->calculate_phenotype(this);
                // this line is to test when male != female, the prob should be 0
                thrust::fill(phenotype[2].begin(), phenotype[2].begin() + size, 1.0f);
                thrust::device_vector<float> zeros(size);
                thrust::fill(zeros.begin(), zeros.end(), 0.0f);
                thrust::transform_if(zeros.begin(), zeros.begin() + size, 
                                  sex.begin(), phenotype[2].begin(), thrust::identity(), isOdd());
                delete genPhenMap;
            }
        }
        ~TestSpeciesAssortNoChance()
        {
            delete neonate;
            delete parents;
        }
};

TEST_F(TestSpeciesAssortNoChance, AssortativeNeonatePairsDemeNoChance) {
    neonate->inherit_genotypes_by_pair(
        parents->parental_pair_probability, 
        parents->males_list, 
        parents->females_list, 
        parents->pair_demes
    );
    // test if kids are from same demes as their moms
    thrust::device_vector<int> mothers_chosen = neonate->get_mothers_chosen();
    ASSERT_NE(neonate->get_mothers_chosen().size(), 0) << "no mothers chosen, which is an unexpected behavior.";
    ASSERT_NE(neonate->kids_deme.size(), 0) << "kids_deme was so poorly constructed and handled.";
    ASSERT_EQ(neonate->get_mothers_chosen().size(), neonate->kids_deme.size()) << "Number of chosen moms must be equal to number of kids.";
    thrust::device_vector<int> mothers_chosen_demes_dev(mothers_chosen.size());
    thrust::gather(mothers_chosen.begin(), mothers_chosen.end(), deme.begin(), mothers_chosen_demes_dev.begin());
    thrust::host_vector<int> kids_deme = neonate->kids_deme;
    thrust::host_vector<int> mothers_chosen_demes = mothers_chosen_demes_dev;
    EXPECT_THAT(mothers_chosen_demes, ::testing::ContainerEq(kids_deme)) << ", in terms of the required equality between mothers' demes and kids' demes.";
    // test if kids are from same demes as their dads
    thrust::device_vector<int> fathers_chosen = neonate->get_fathers_chosen();
    ASSERT_NE(neonate->get_fathers_chosen().size(), 0) << "no fathers chosen, which is an unexpected behavior.";
    ASSERT_EQ(neonate->get_fathers_chosen().size(), neonate->kids_deme.size()) << "Number of chosen dads must be equal to number of kids.";
    thrust::device_vector<int> fathers_chosen_demes_dev(fathers_chosen.size());
    thrust::gather(fathers_chosen.begin(), fathers_chosen.end(), deme.begin(), fathers_chosen_demes_dev.begin());
    thrust::host_vector<int> fathers_chosen_demes = fathers_chosen_demes_dev;
    EXPECT_THAT(fathers_chosen_demes, ::testing::ContainerEq(kids_deme)) << ", in terms of the required equality between fathers' demes and kids' demes.";
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

        void SetUp() override 
        {
            for (int i = 0; i < demeParameters->check_number_of_demes(); i++)
            {
                ASSERT_NE(deme_sizes[i], 0) << "Index " << i << " is 0";
            }
        }

        void initializeDemes()
        {
            int demeSizeArr[Num_Demes];
            distributeIndvsIntoDemes(demeSizeArr, size, Num_Demes);
            thrust::fill(max_deme_sizes.begin(), max_deme_sizes.begin() + Num_Demes, maxsize/Num_Demes);
            int sum_so_far = 0;
            for (int i = 0; i < Num_Demes; i++)
            {
                int start = sum_so_far;
                int end = start + demeSizeArr[i];
                thrust::fill(deme.begin() + start, deme.begin() + end, i);
                sum_so_far += demeSizeArr[i];
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
    thrust::device_vector<int> maternal_id_chosen(neonate->Total_Number_of_Neonates);
    thrust::copy(maternal_id.begin() + neonate->previous_pop_size, 
        maternal_id.begin() + neonate->previous_pop_size + neonate->Total_Number_of_Neonates, 
        maternal_id_chosen.begin());
    thrust::device_vector<int> mothers_deme(neonate->Total_Number_of_Neonates);
    thrust::gather(maternal_id_chosen.begin(), maternal_id_chosen.end(), deme.begin(), mothers_deme.begin());
    thrust::host_vector<int> mothers_deme_copy = mothers_deme;
    thrust::host_vector<int> kids_deme_copy = neonate->kids_deme;
    EXPECT_THAT(mothers_deme_copy, ::testing::ContainerEq(kids_deme_copy));
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
