#include "parser.hpp"
#include "tensor.hpp"
#include "neural.hpp"

#define BATCH_SIZE 20

int main (int argc, char** argv) {
    std::unordered_map<std::string, std::vector<double>> data;
    data = parse_csv<double>("diabetes.csv");
    print_map(data);



    Network* net = new Network();
    net->add(new Dense(8, 4, BATCH_SIZE));
    net->add(new Activation(4, BATCH_SIZE, act::tanh, der::tanh));

    net->add(new Dense(4, 1, BATCH_SIZE));
    net->add(new Activation(1, BATCH_SIZE, act::sigmoid, der::sigmoid));


    Model mod(net, err::binary);
    mod.params(400,40,BATCH_SIZE,0.1,0.8);
    mod.train(data, {"Outcome"});
    return 0;
}
