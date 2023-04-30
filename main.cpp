#include "parser.hpp"
#include "tensor.hpp"
#include "neural.hpp"
#include "progress.hpp"

#define BATCH_SIZE 120

int main (int argc, char** argv) {
    signal(SIGINT, restore_terminal);
    Dictionary<std::vector<double>> data = parse_general_csv<double>("datasets/IRIS.csv");
    make_one_hot(data, "species", {"setosa", "versicolor", "virginica"});
    print_dictionary(data);

    Model mod(err::categorical, BATCH_SIZE, 4);
    mod.add(LAYER::DENSE, ACTIVATION::TANH, 4);
    mod.add(LAYER::DENSE, ACTIVATION::SOFTMAX, 3);
    mod.params(4000,1,0.01,0.8);
    mod.train(data, {"setosa", "versicolor", "virginica"});

    return 0;
}
