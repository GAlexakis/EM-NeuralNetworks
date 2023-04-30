#include "parser.hpp"
#include "tensor.hpp"
#include "neural.hpp"
#include "progress.hpp"

#define BATCH_SIZE 100

int main (int argc, char** argv) {
    signal(SIGINT, restore_terminal);
    Dictionary<std::vector<double>> data = parse_csv<double>("datasets/diabetes.csv");

    Model mod(err::regression, BATCH_SIZE, 8);
    mod.add(LAYER::DENSE, ACTIVATION::TANH, 6);
    mod.add(LAYER::DENSE, ACTIVATION::TANH, 4);
    mod.add(LAYER::DENSE, ACTIVATION::TANH, 2);
    mod.add(LAYER::DENSE, ACTIVATION::SIGMOID, 1);
    mod.params(100,80,0.001,0.8);
    mod.train(data, {"Outcome"});

    return 0;
}
