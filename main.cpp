#include "parser.hpp"
#include "tensor.hpp"
#include "neural.hpp"

void fun(Tensor<double>* t1) {
    *t1 = exp(*t1);
}

int main (int argc, char** argv) {
    std::vector<size_t> dims = {1,5};
    std::vector<double> vals = {1,2,3,4,5};
    Tensor<double>* t1 = new Tensor(vals, dims);
    Network* net = new Network();
    net->add(new Activation(5, 1, fun, fun));
    t1->print();
    net->forward(t1);
    t1->print();
    net->backwards(t1);
    t1->print();
    net->fix(0.5);
    delete t1;
    delete net;
    return 0;
}
