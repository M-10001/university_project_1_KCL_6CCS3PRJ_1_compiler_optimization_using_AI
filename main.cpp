#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <unordered_map>

// Plotting.
#include "matplotlibcpp.h"

using namespace std;

// Plotting.
namespace plt = matplotlibcpp;

// Location of compiler optimization flags information file.
const string FLAG_FILE = "compiler_optimization_flags_information.csv";

// Population size of chromosomes created. Must not be lower than 2.
const int POPULATION_SIZE = 20;

// Used to determine off value for Gene in Chromosome.
const int OFF_GENE_VALUE = -1;

// Number of passes that will not be used in run time calculation.
const int UNUSED_RUN_PASSES = 0;

// Used for the number of times the output file will be ran. Must not be lower than or equal to UNUSED_RUN_PASSES.
const int RUN_PASSES = 5;

// Probability of a gene variant to mutate to a different variant.
const double MUTATION_PROBABILITY = 0.1;

// Used for increasing or decreasing weight of a Gene variant.
const int VARIANT_WEIGHT_INCREMENT = 1;

// Seed for randomizer.
const int RANDOMIZER_SEED = 0;

// The maximum number of passes in the main loop.
const int MAX_PASSES = 100;

// Used to check against for updating fixed variants based on number of passes occurred.
const int FIXED_VARIANTS_UPDATE_PASS = 20;

// Flags to compare against for current algorithm's improvement.
const vector<string> FLAGS_TO_COMPARE_AGAINST = {" -O1", " -O2", " -O3"};

// The number of passes the comparison flags will be recompiled and ran over.
const int COMPARISON_FLAGS_PASSES = 3;

// Value needed in calculation to determine the threshold if a chromosome run time is an outlier.
const double A = 1.1;

// Value needed in calculation to determine the threshold if a chromosome run time is an outlier.
const double B = 15.0;

// Value needed in calculation to determine the threshold if a chromosome run time is an outlier.
const double C = 1.0;

// Value in fitness calculation. Must be between 0 and 1.
const double D = 0.99;

// Value in selecting Genes variants to be fixed on.
const int E = 2;

struct GeneVariant {
    string flag;
    int weight = 0;

    GeneVariant(string flag) : flag(flag){}
};

struct Gene {
    vector<GeneVariant> variants;

    Gene(vector<GeneVariant> variants) : variants(variants){}
};

struct Chromosome {
    vector<int> on_off_for_variants;
    double run_time = 0.0;
    double fitness = 0.0;

    Chromosome(vector<int> on_off_for_variants) : on_off_for_variants(on_off_for_variants){}
};

vector<string> get_compilation_command_parts(){
    vector<string> command_parts;
    string command;
    cout << "Input compilation command:" << endl;
    getline(cin, command);
    int split_position = command.find(" ");
    command_parts.push_back(command.substr(0, split_position));
    command_parts.push_back(command.substr(split_position, command.size() - split_position));
    return command_parts;
}

string get_output_file_path_with_inputs(){
    cout << "Input output file path with inputs:" << endl;
    string output_file_path_with_inputs;
    getline(cin, output_file_path_with_inputs);
    return output_file_path_with_inputs;
}

string get_output_file_path(const string & output_file_path_with_inputs){
    int i = 0;

    while (output_file_path_with_inputs[i] == ' '){
        i = i + 1;
    }

    while (output_file_path_with_inputs[i] != ' '){
        i = i + 1;
    }

    return output_file_path_with_inputs.substr(0, i);
}

void combine_compilation_command_parts_to_flags(string & flags, const vector<string> & compilation_command_parts){
    flags = compilation_command_parts[0] + flags + compilation_command_parts[1];
}

void compile_with_command(const string & command){
    FILE * pipe = popen(command.c_str(), "r");
    int compile_status = pclose(pipe);

    if (compile_status != 0){
        throw runtime_error("\nError when compiling with command: " + command);
    }
}

double get_run_time(const string & file_path_with_inputs, const string & compilation_command){
    auto start_time = chrono::high_resolution_clock::now();
    FILE * pipe = popen(file_path_with_inputs.c_str(), "r");
    double run_status = pclose(pipe);
    auto end_time = chrono::high_resolution_clock::now();

    if (run_status != 0){
        throw runtime_error(
            "\nError when running file and any added inputs: " +
            file_path_with_inputs +
            "\nCreated by command: " +
            compilation_command
        );
    }

    chrono::duration<double, milli> time_duration = end_time - start_time;
    return time_duration.count();
}

double compile_and_get_avg_run_time(const string & compilation_command, const string & output_file_path_with_inputs){
    compile_with_command(compilation_command);
    double total_run_time = 0.0;
    int i = 0;

    while (i < RUN_PASSES){
        double run_time = get_run_time(output_file_path_with_inputs, compilation_command);

        if (i >= UNUSED_RUN_PASSES){
            total_run_time = total_run_time + run_time;
        }

        i++;
    }

    remove(get_output_file_path(output_file_path_with_inputs).c_str());
    return (total_run_time / (double) (i - UNUSED_RUN_PASSES));
}

void set_genes(vector<Gene> & genes){
    ifstream file_reader(FLAG_FILE);
    string line = "";

    while (getline(file_reader, line)){
        stringstream line_stream(line);
        vector<GeneVariant> variants;

        if (line.find(',') == string::npos){
            variants.push_back(GeneVariant(" " + line));
        } else {
            string start_word = "";
            getline(line_stream, start_word, ',');
            string word = "";

            while (getline(line_stream, word, ',')){
                variants.push_back(GeneVariant(" " + start_word + word));
            }
        }

        genes.push_back(Gene(variants));
    }

    file_reader.close();
}

string get_chromosome_compilation_command(
    const Chromosome & chromosome,
    const vector<Gene> & genes,
    const vector<string> & compliation_command_parts
){
    string compilation_command = "";

    for (int i = 0; i < genes.size(); i++){
        if (chromosome.on_off_for_variants[i] != OFF_GENE_VALUE){
            compilation_command = compilation_command + genes[i].variants[chromosome.on_off_for_variants[i]].flag;
        }
    }

    combine_compilation_command_parts_to_flags(compilation_command, compliation_command_parts);
    return compilation_command;
}

void set_run_time(
    Chromosome & chromosome,
    const vector<Gene> & genes,
    const vector<string> & compliation_command_parts,
    const string & output_file_path_with_inputs
){
    string compilation_command = get_chromosome_compilation_command(chromosome, genes, compliation_command_parts);
    chromosome.run_time = compile_and_get_avg_run_time(compilation_command, output_file_path_with_inputs);
}

void set_run_times(
    vector<Chromosome> & chromosomes,
    const vector<Gene> & genes,
    const vector<string> & compliation_command_parts,
    const string & output_file_path_with_inputs
){
    for (int i = 0; i < chromosomes.size(); i++){
        set_run_time(
            chromosomes[i],
            genes,
            compliation_command_parts,
            output_file_path_with_inputs
        );
    }
}

double get_chromosomes_run_time_avg(const vector<Chromosome> & chromosomes){
    double run_time_sum = 0.0;
    int i = 0;

    while (i < chromosomes.size()){
        run_time_sum = run_time_sum + chromosomes[i].run_time;
        i++;
    }

    return (run_time_sum / (double) i);
}

double get_chromosomes_run_time_s_d(const vector<Chromosome> & chromosomes, const double & mean){
    double standard_deviation = 0.0;

    for (int i = 0; i < chromosomes.size(); i++){
        standard_deviation = standard_deviation + pow(chromosomes[i].run_time - mean, 2);
    }

    standard_deviation = sqrt(standard_deviation / (double) chromosomes.size());

    return standard_deviation;
}

void set_run_times_for_outliers(
    vector<Chromosome> & chromosomes,
    const vector<Gene> genes,
    const vector<string> & compliation_command_parts,
    const string & output_file_path_with_inputs
){
    double mean = get_chromosomes_run_time_avg(chromosomes);
    double standard_deviation = get_chromosomes_run_time_s_d(chromosomes, mean);
    double threshold = (B * (1 - pow(A, -standard_deviation))) + C;

    int j = 0;

    for (int i = 0; i < chromosomes.size(); i++){
        if (abs(chromosomes[i].run_time - mean) > threshold){
            set_run_time(
                chromosomes[i],
                genes,
                compliation_command_parts,
                output_file_path_with_inputs
            );

            j++;
        }
    }

    cout << "Number of re-run chromosomes: " << j << endl;
}

Chromosome get_min_run_time_chromosome(const vector<Chromosome> & chromosomes){
    Chromosome min_run_time_chromosome(chromosomes[0]);

    for (int i = 1; i < chromosomes.size(); i++){
        if (chromosomes[i].run_time < min_run_time_chromosome.run_time){
            min_run_time_chromosome = chromosomes[i];
        }
    }

    return min_run_time_chromosome;
}

void set_min_run_time_chromosome(Chromosome & min_run_time_chromosome, const vector<Chromosome> & chromosomes){
    for (int i = 0; i < chromosomes.size(); i++){
        if (chromosomes[i].run_time < min_run_time_chromosome.run_time){
            min_run_time_chromosome = chromosomes[i];
        }
    }
}

void set_fitnesses(vector<Chromosome> & chromosomes, const int & min_run_time){
    for (int i = 0; i < chromosomes.size(); i++){
        chromosomes[i].fitness = 1 / ((chromosomes[i].run_time - (min_run_time * D)) / (min_run_time * D));
    }
}

void set_initial_population(
    vector<Chromosome> & population,
    const vector<Gene> & genes
){
    for (int i = 0; i < POPULATION_SIZE; i++){
        vector<int> on_off_for_variants;

        for (int j = 0; j < genes.size(); j++){
            int on_off_variant_randomizer = rand() % (genes[j].variants.size() + 1);

            if (on_off_variant_randomizer < genes[j].variants.size()){
                on_off_for_variants.push_back(on_off_variant_randomizer);
            } else {
                on_off_for_variants.push_back(OFF_GENE_VALUE);
            }
        }

        Chromosome chromosome(on_off_for_variants);
        population.push_back(chromosome);
    }
}

void choose_c_1_and_c_2(const Chromosome * & c_1, const Chromosome * & c_2, const vector<Chromosome> & chromosomes){
    double fitness_sum = 0.0;

    for (int i = 0; i < chromosomes.size(); i++){
        fitness_sum = fitness_sum + chromosomes[i].fitness;
    }

    double chromosome_randomizer = rand() / (double) RAND_MAX;
    double probability_accumulator = 0.0;
    int i = 0;

    while ((i < chromosomes.size()) && (c_1 == nullptr)){
        double chromosome_probability = chromosomes[i].fitness / fitness_sum;
        probability_accumulator = probability_accumulator + chromosome_probability;

        if ((chromosome_probability > 0.0) && (chromosome_randomizer <= probability_accumulator)){
            c_1 = &chromosomes[i];
        }

        i++;
    }

    c_2 = c_1;

    while (c_2 == c_1){
        chromosome_randomizer = rand() / (double) RAND_MAX;
        probability_accumulator = 0.0;
        c_2 = nullptr;
        i = 0;

        while ((i < chromosomes.size()) && (c_2 == nullptr)){
            double chromosome_probability = chromosomes[i].fitness / fitness_sum;
            probability_accumulator = probability_accumulator + chromosome_probability;

            if ((chromosome_probability > 0.0) && (chromosome_randomizer <= probability_accumulator)){
                c_2 = &chromosomes[i];
            }

            i++;
        }
    }
}

void recombine_r_c_1_and_r_c_2(Chromosome & r_c_1, Chromosome & r_c_2){
    int locus = rand() % r_c_1.on_off_for_variants.size();

    if ((rand() % 2) == 0){
        for (int i = 0; i < locus; i++){
            int temp = r_c_1.on_off_for_variants[i];
            r_c_1.on_off_for_variants[i] = r_c_2.on_off_for_variants[i];
            r_c_2.on_off_for_variants[i] = temp;
        }
    } else {
        for (int i = locus; i < r_c_1.on_off_for_variants.size(); i++){
            int temp = r_c_1.on_off_for_variants[i];
            r_c_1.on_off_for_variants[i] = r_c_2.on_off_for_variants[i];
            r_c_2.on_off_for_variants[i] = temp;
        }
    }
}

void mutate_gene(Chromosome & chromosome, const int & gene_pos, const Gene & gene){
    int min_weight = 0;
    int weight_sum = 0;
    int i = 0;

    while (i < gene.variants.size()){
        weight_sum = weight_sum + gene.variants[i].weight;

        if (gene.variants[i].weight < min_weight){
            min_weight = gene.variants[i].weight;
        }

        i++;
    }

    double adder = 1;

    if (min_weight < 0){
        adder = adder - min_weight;
    }

    double normalized_weight_sum = (adder * (i + 1)) + weight_sum;
    double mutation_randomizer = rand() / (double) RAND_MAX;
    double probability_accumulator = 0.0;
    bool mutated = false;
    i = 0;

    while ((i < gene.variants.size()) && (mutated == false)){
        double mutation_probability = (adder + gene.variants[i].weight) / normalized_weight_sum;
        probability_accumulator = probability_accumulator + mutation_probability;

        if ((mutation_probability > 0.0) && (mutation_randomizer <= probability_accumulator)){
            chromosome.on_off_for_variants[gene_pos] = i;
            mutated = true;
        }

        i++;
    }

    if (mutated == false){
        chromosome.on_off_for_variants[gene_pos] = OFF_GENE_VALUE;
    }
}

void mutate_chromosome(Chromosome & chromosome, const vector<Gene> & genes){
    for (int i = 0; i < genes.size(); i++){
        double mutation_randomizer = rand() / (double) RAND_MAX;

        if ((MUTATION_PROBABILITY > 0.0) && (mutation_randomizer <= MUTATION_PROBABILITY)){
            mutate_gene(chromosome, i, genes[i]);
        }
    }
}

void fix_chromosome(Chromosome & chromosome, const unordered_map<int, int> & fixed_variants_pos){
    for (auto const& i : fixed_variants_pos){
        chromosome.on_off_for_variants[i.first] = i.second;
    }
}

void remove_chromosomes(vector<Chromosome> & chromosomes, const int & removal_amount){
    double reverse_fitness_sum = 0;

    for (int i = 0; i < chromosomes.size(); i++){
        reverse_fitness_sum = reverse_fitness_sum + (1 / chromosomes[i].fitness);
    }

    for (int i = 0; i < removal_amount; i++){
        double chromosome_randomizer = rand() / (double) RAND_MAX;
        double probability_accumulator = 0.0;
        bool removed = false;
        int j = 0;

        while ((removed == false) && (j < chromosomes.size())){
            double reverse_chromosome_probability = (1 / chromosomes[j].fitness) / reverse_fitness_sum;
            probability_accumulator = probability_accumulator + reverse_chromosome_probability;

            if ((reverse_chromosome_probability > 0.0) && (chromosome_randomizer <= probability_accumulator)){
                reverse_fitness_sum = reverse_fitness_sum - (1 / chromosomes[j].fitness);
                chromosomes.erase(chromosomes.begin() + j);
                removed = true;
            }

            j++;
        }
    }
}

void increment_variant_weights(
    Gene & gene,
    const int & variant_pos,
    const int & comparison_variant_pos,
    const bool & improved
){
    if (improved == true){
        if ((variant_pos != OFF_GENE_VALUE) && (comparison_variant_pos == OFF_GENE_VALUE)){
            gene.variants[variant_pos].weight = gene.variants[variant_pos].weight + VARIANT_WEIGHT_INCREMENT;
        } else if ((variant_pos == OFF_GENE_VALUE) && (comparison_variant_pos != OFF_GENE_VALUE)){
            gene.variants[comparison_variant_pos].weight = gene.variants[comparison_variant_pos].weight - VARIANT_WEIGHT_INCREMENT;
        } else if ((variant_pos != OFF_GENE_VALUE) && (comparison_variant_pos != OFF_GENE_VALUE) && (variant_pos != comparison_variant_pos)) {
            gene.variants[variant_pos].weight = gene.variants[variant_pos].weight + VARIANT_WEIGHT_INCREMENT;
            gene.variants[comparison_variant_pos].weight = gene.variants[comparison_variant_pos].weight - VARIANT_WEIGHT_INCREMENT;
        }
    } else {
        if ((variant_pos != OFF_GENE_VALUE) && (comparison_variant_pos == OFF_GENE_VALUE)){
            gene.variants[variant_pos].weight = gene.variants[variant_pos].weight - VARIANT_WEIGHT_INCREMENT;
        } else if ((variant_pos == OFF_GENE_VALUE) && (comparison_variant_pos != OFF_GENE_VALUE)){
            gene.variants[comparison_variant_pos].weight = gene.variants[comparison_variant_pos].weight + VARIANT_WEIGHT_INCREMENT;
        } else if ((variant_pos != OFF_GENE_VALUE) && (comparison_variant_pos != OFF_GENE_VALUE) && (variant_pos != comparison_variant_pos)) {
            gene.variants[variant_pos].weight = gene.variants[variant_pos].weight - VARIANT_WEIGHT_INCREMENT;
            gene.variants[comparison_variant_pos].weight = gene.variants[comparison_variant_pos].weight + VARIANT_WEIGHT_INCREMENT;
        }
    }
}

void increment_variants_weights_and_add_candidate(
    vector<Gene> & genes,
    vector<Chromosome> & candidate_list,
    const Chromosome & chromosome,
    const Chromosome & comparison_chromosome
){
    bool improved = false;

    if (chromosome.fitness > comparison_chromosome.fitness){
        candidate_list.push_back(chromosome);
        improved = true;
    }

    for (int i = 0; i < genes.size(); i++){
        increment_variant_weights(genes[i], chromosome.on_off_for_variants[i], comparison_chromosome.on_off_for_variants[i], improved);
    }
}

void set_fixed_variants_pos(
    unordered_map<int, int> & fixed_variants_pos,
    const vector<Chromosome> & candidate_list
){
    unordered_map<int, unordered_map<int, int>> variant_pos_tracker;

    for (int i = 0; i < candidate_list.size(); i++){
        for (int j = 0; j < candidate_list[i].on_off_for_variants.size(); j++){
            int variant_pos = candidate_list[i].on_off_for_variants[j];

            if (
                (variant_pos_tracker.find(j) != variant_pos_tracker.end()) &&
                (variant_pos_tracker[j].find(variant_pos) != variant_pos_tracker[j].end())
            ){
                variant_pos_tracker[j][variant_pos] = variant_pos_tracker[j][variant_pos] + 1;
            }
            else if (variant_pos != OFF_GENE_VALUE) {
                unordered_map<int, int> variant_occurrence_tracker = {{variant_pos, 1}};
                variant_pos_tracker[j] = variant_occurrence_tracker;
            }
        }
    }

    fixed_variants_pos.clear();
    int fixing_requirement = candidate_list.size() - E;

    for (auto const& i : variant_pos_tracker) {
        for (auto const& j : i.second) {
            if (
                (j.second >= fixing_requirement) &&
                (
                    (fixed_variants_pos.find(i.first) == fixed_variants_pos.end()) ||
                    (variant_pos_tracker[i.first][fixed_variants_pos[i.first]] < j.second)
                )
            ){
                fixed_variants_pos[i.first] = j.first;
            }
        }
    }
}

void set_comparison_flag_on_min_run_time(
    string & flag,
    double & min_run_time,
    const vector<string> & compilation_command_parts,
    const string & output_file_name
){
    vector<double> run_times_sum;

    for (int i = 0; i < FLAGS_TO_COMPARE_AGAINST.size(); i++){
        run_times_sum.push_back(0.0);

        for (int j = 0; j < COMPARISON_FLAGS_PASSES; j++){
            string compilation_command = FLAGS_TO_COMPARE_AGAINST[i];
            combine_compilation_command_parts_to_flags(compilation_command, compilation_command_parts);
            run_times_sum[i] = run_times_sum[i] + compile_and_get_avg_run_time(compilation_command, output_file_name);
        }
    }

    for (int i = 0; i < FLAGS_TO_COMPARE_AGAINST.size(); i++){
        double run_time = run_times_sum[i] / COMPARISON_FLAGS_PASSES;

        if ((i == 0) || (run_time < min_run_time)){
            min_run_time = run_time;
            flag = FLAGS_TO_COMPARE_AGAINST[i].substr(1 , FLAGS_TO_COMPARE_AGAINST[i].size() - 1);
        }
    }
}

// Plotting.
void plot_output(
    const vector<double> & run_times_avg,
    const vector<double> & run_times_s_d,
    const double & comparison_flag_run_time
){
    vector<int> x_axis;

    for (int i = 0; i < run_times_avg.size(); i++){
        x_axis.push_back(i);
    }

    plt::figure();
    plt::subplot(1, 3, 1);
    plt::plot(x_axis, run_times_avg, "r--");
    plt::title("Run time average of population:");
    plt::xlabel("Pass numbers");
    plt::ylabel("Run times average (ms)");

    plt::subplot(1, 3, 2);
    plt::plot(x_axis, run_times_s_d, "r--");
    plt::title("Run time standard deviation of population:");
    plt::xlabel("Pass numbers");
    plt::ylabel("Run times standard deviation (ms)");

    vector<double> avg_run_times_percentage_increase_from_comparison;

    for (int i = 0; i < run_times_avg.size(); i++){
        avg_run_times_percentage_increase_from_comparison.push_back(100 * ((comparison_flag_run_time - run_times_avg[i]) / comparison_flag_run_time));
    }

    plt::subplot(1, 3, 3);
    plt::plot(x_axis, avg_run_times_percentage_increase_from_comparison, "r--");
    plt::title("Average run time improvement from comparison flags:");
    plt::xlabel("Pass numbers");
    plt::ylabel("Percentage (%)");
    plt::show();
}

int main(){
    vector<string> compilation_command_parts = get_compilation_command_parts();

    cout << "------" << endl;

    string output_file_path_with_inputs = get_output_file_path_with_inputs();

    cout << "------" << endl;

    auto start_time = chrono::high_resolution_clock::now();

    cout << "Time started at 0ms" << endl;

    cout << "------" << endl;
    cout << "Initialization started" << endl;

    srand(RANDOMIZER_SEED);

    vector<Gene> genes;
    set_genes(genes);

    vector<Chromosome> population;

    set_initial_population(
        population,
        genes
    );

    set_run_times(
        population,
        genes,
        compilation_command_parts,
        output_file_path_with_inputs
    );

    set_run_times_for_outliers(
        population,
        genes,
        compilation_command_parts,
        output_file_path_with_inputs
    );

    vector<double> run_times_avg = {get_chromosomes_run_time_avg(population)};
    vector<double> run_times_s_d = {get_chromosomes_run_time_s_d(population, run_times_avg.back())};

    Chromosome min_run_time_chromosome(get_min_run_time_chromosome(population));
    set_fitnesses(population, min_run_time_chromosome.run_time);

    cout << "Current run time average: " << run_times_avg.back() << "ms" << endl;
    cout << "Current run time standard deviation: " << run_times_s_d.back() << "ms" << endl;

    vector<Chromosome> candidate_list;
    unordered_map<int, int> fixed_variants_pos;
    int since_significant_update_pass = 0;
    int pass_number = 0;

    cout << "------" << endl;
    cout << "Main loop started" << endl;

    while (pass_number < MAX_PASSES){
        cout << "------" << endl;

        const Chromosome * c_1 = nullptr;
        const Chromosome * c_2 = nullptr;
        choose_c_1_and_c_2(c_1, c_2, population);

        Chromosome r_c_1(*c_1);
        Chromosome r_c_2(*c_2);
        recombine_r_c_1_and_r_c_2(r_c_1, r_c_2);

        Chromosome m_c_1(r_c_1);
        Chromosome m_c_2(r_c_2);
        mutate_chromosome(m_c_1, genes);
        mutate_chromosome(m_c_2, genes);

        Chromosome g_c_1(r_c_1);
        Chromosome g_c_2(r_c_2);
        fix_chromosome(g_c_1, fixed_variants_pos);
        fix_chromosome(g_c_2, fixed_variants_pos);

        vector<Chromosome> mutated_chromosomes = {m_c_1, m_c_2, g_c_1, g_c_2, r_c_2, r_c_1};

        set_run_times(
            mutated_chromosomes,
            genes,
            compilation_command_parts,
            output_file_path_with_inputs
        );

        set_min_run_time_chromosome(min_run_time_chromosome, mutated_chromosomes);
        set_fitnesses(mutated_chromosomes, min_run_time_chromosome.run_time);

        r_c_1 = mutated_chromosomes.back();
        mutated_chromosomes.pop_back();
        r_c_2 = mutated_chromosomes.back();
        mutated_chromosomes.pop_back();

        increment_variants_weights_and_add_candidate(genes, candidate_list, mutated_chromosomes[0], r_c_1);
        increment_variants_weights_and_add_candidate(genes, candidate_list, mutated_chromosomes[1], r_c_2);
        increment_variants_weights_and_add_candidate(genes, candidate_list, mutated_chromosomes[2], r_c_1);
        increment_variants_weights_and_add_candidate(genes, candidate_list, mutated_chromosomes[3], r_c_2);

        population.insert(
            population.end(),
            mutated_chromosomes.begin(),
            mutated_chromosomes.end()
        );

        set_run_times_for_outliers(
            population,
            genes,
            compilation_command_parts,
            output_file_path_with_inputs
        );

        min_run_time_chromosome = get_min_run_time_chromosome(population);
        set_fitnesses(population, min_run_time_chromosome.run_time);

        remove_chromosomes(population, (population.size() - POPULATION_SIZE));

        run_times_avg.push_back(get_chromosomes_run_time_avg(population));
        run_times_s_d.push_back(get_chromosomes_run_time_s_d(population, run_times_avg.back()));

        min_run_time_chromosome = get_min_run_time_chromosome(population);
        set_fitnesses(population, min_run_time_chromosome.run_time);

        if (((pass_number + 1) % FIXED_VARIANTS_UPDATE_PASS) == 0){
            set_fixed_variants_pos(fixed_variants_pos, candidate_list);
            candidate_list.clear();
        }

        cout << "Current run time average: " << run_times_avg.back() << "ms" << endl;
        cout << "Current run time standard deviation: " << run_times_s_d.back() << "ms" << endl;
        cout << "Completed pass: " << pass_number << endl;

        pass_number++;
    }

    string compilation_command = "";
    combine_compilation_command_parts_to_flags(compilation_command, compilation_command_parts);
    double base_run_time = compile_and_get_avg_run_time(compilation_command, output_file_path_with_inputs);

    string comparison_flag = "";
    double comparison_flag_run_time = 0.0;

    set_comparison_flag_on_min_run_time(
        comparison_flag,
        comparison_flag_run_time,
        compilation_command_parts,
        output_file_path_with_inputs
    );

    compilation_command = get_chromosome_compilation_command(min_run_time_chromosome, genes, compilation_command_parts);
    double min_run_time_chromosome_re_run_time = compile_and_get_avg_run_time(compilation_command, output_file_path_with_inputs);

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> time_duration = end_time - start_time;
    double total_time_taken = time_duration.count();
    string time_classifier = "ms";

    if (total_time_taken > 3600000){
        total_time_taken = total_time_taken / 3600000;
        time_classifier = "hr";
    } else if (total_time_taken > 60000){
        total_time_taken = total_time_taken / 60000;
        time_classifier = "min";
    } else if (total_time_taken > 1000){
        total_time_taken = total_time_taken / 1000;
        time_classifier = "s";
    }

    cout << "------" << endl;
    cout << "Final outputs:" << endl;

    cout << "Base run time: " << base_run_time << "ms" << endl;

    cout << "Comparison flag: " << comparison_flag << " and run time: " << comparison_flag_run_time << "ms" << endl;

    cout << "Min run time out of current population: " << min_run_time_chromosome.run_time << "ms" << endl;
    cout << "Min run time chromosome out of current population compilation command:" << endl;
    cout << compilation_command << endl;

    cout << "Re-run time of min run time chromosome: " << min_run_time_chromosome_re_run_time << "ms" << endl;

    cout << "Total time taken: " << total_time_taken << time_classifier << endl;

    // Plotting.
    plot_output(run_times_avg, run_times_s_d, comparison_flag_run_time);

    cout << "------" << endl;
    cout << "Press enter to continue" << endl;
    cin.ignore();
    return 0;
}
