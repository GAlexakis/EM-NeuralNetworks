#include "parser.hpp"
#include "tensor.hpp"
#include "neural.hpp"

void fun(Tensor<double>* t1) {
    *t1 = exp(*t1);
}

int main (int argc, char** argv) {
    std::vector<size_t> dims = {1,5};
    std::vector<double> vals1 = {-2,-1,0,1,2};
    std::vector<double> vals2 = {0,1,0,1,0};
    Tensor<double>* t1 = new Tensor(vals1, dims);
    Tensor<double> t2 = Tensor(vals2, dims);
    // Network* net = new Network();
    // net->add(new Activation(5, 2, act::softmax, fun));
    t1->print();

    (*t1 > 0).print();
    der::sigmoid(t1);
    // net->forward(t1);
    t1->print();
    // std::cout << "MIN " << std::numeric_limits<double>::min() << '\n';
    delete t1;
    // delete net;
    return 0;
}
