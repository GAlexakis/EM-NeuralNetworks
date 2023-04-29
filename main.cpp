#include "parser.hpp"
#include "tensor.hpp"
#include "neural.hpp"

#define BATCH_SIZE 100

int main (int argc, char** argv) {
    std::unordered_map<std::string, std::vector<double>> data;
    // data = parse_csv<double>("diabetes.csv");
    data.insert({"value1", {}});
    data.insert({"value2", {}});
    data.insert({"Outcome", {}});
    for (double i = 0; i < 400; i++) {

        double val1 = 8*rand()/(double)(RAND_MAX) - 4;
        double val2 = 8*rand()/(double)(RAND_MAX) - 4;

        double res = 0.7*exp(cos(M_PIl*val1)) + 0.3*cos(2*M_PIl*val2);

        data["value1"].push_back(val1);
        data["value2"].push_back(val2);
        data["Outcome"].push_back(res);
    }



    Network* net = new Network();

    net->add(new Dense(2, 4, BATCH_SIZE));
    net->add(new Activation(4, BATCH_SIZE, act::tanh, der::tanh));

    net->add(new Dense(4, 8, BATCH_SIZE));
    net->add(new Activation(8, BATCH_SIZE, act::tanh, der::tanh));

    net->add(new Dense(8, 4, BATCH_SIZE));
    net->add(new Activation(4, BATCH_SIZE, act::tanh, der::tanh));

    net->add(new Dense(4, 2, BATCH_SIZE));
    net->add(new Activation(2, BATCH_SIZE, act::tanh, der::tanh));

    net->add(new Dense(2, 1, BATCH_SIZE));
    net->add(new Activation(1, BATCH_SIZE, act::sigmoid, der::sigmoid));


    Model mod(net, err::regression);
    mod.params(2000,1,BATCH_SIZE,0.001,0.8);
    mod.train(data, {"Outcome"});


    double user_input;
    while (true) {
        std::cin >> user_input;
        mod.predict({user_input});
    }
    return 0;
}
