#include "parser.hpp"
#include "tensor.hpp"
#include "neural.hpp"

Dense::Dense (size_t input_size, size_t output_size, size_t batch_size) {
    weights = new Tensor<double>(input_size, output_size);
    biases  = new Tensor<double>(1, output_size);
    inputs  = new Tensor<double>(batch_size, input_size);
    errors  = new Tensor<double>(batch_size, output_size);

    for( size_t i = 0; i < input_size; i++) {
            for( size_t j = 0; j < output_size; j++) {
                (*weights)(i, j) = rand()/(double)(RAND_MAX) - 0.5;
            }
        }
        for( size_t i = 0; i < output_size; i++) {
            (*biases)(0, i) = 0;
        }
}
Dense::~Dense () {
    delete weights;
    delete biases;
    delete inputs;
    delete errors;
}
void Dense::forward (Tensor<double>* propagator) {
    *inputs = *propagator;
    *propagator = mul(*propagator, *weights) + *biases;
}
void Dense::backwards (Tensor<double>* propagator) {
    *errors = *propagator;
    *propagator = mul(*propagator, ~(*weights));
}
void Dense::predict (Tensor<double>* propagator) {
    *propagator = mul(*propagator, *weights) + *biases;
}
void Dense::fix (double learning_rate) {
    weights->print();
    biases->print();
    *weights = *weights - mul(~(*inputs), *errors)*learning_rate/(double)errors->dims()[0];
    *biases = *biases - (*errors)[0]*learning_rate/(double)errors->dims()[0];
    weights->print();
    biases->print();
}

Activation::Activation (size_t input_size, size_t batch_size, std::function<void(Tensor<double>*)> func, std::function<void(Tensor<double>*)> der) {
    inputs  = new Tensor<double>(batch_size, input_size);
    this->func = func;
    this->der = der;
}
Activation::~Activation () {
    func = nullptr;
    der = nullptr;
    delete inputs;
}
void Activation::forward (Tensor<double>* propagator) {
    *inputs = *propagator;
    func(propagator);
}
void Activation::backwards (Tensor<double>* propagator) {
    der(inputs);
    *propagator = (*propagator) * (*inputs);
}
void Activation::predict (Tensor<double>* propagator) {
    func(propagator);
}
void Activation::fix (double learning_rate) {

}

Network::Network () {
    layers = std::vector<Layer*>{};
}
Network::~Network () {
    for (auto layer : layers) delete layer;
 }
void Network::add (Layer* l) {
    layers.push_back(l);
}
void Network::forward (Tensor<double>* propagator) {
    for (auto layer : layers) layer->forward(propagator);
}
void Network::backwards (Tensor<double>* propagator) {
    for(size_t i = 0; i < layers.size(); i++)
        layers[layers.size() - i - 1]->backwards(propagator);
}
void Network::predict (Tensor<double>* propagator) {
    for (auto layer : layers) layer->predict(propagator);
}
void Network::fix (double learning_rate) {
    for (auto layer : layers) layer->fix(learning_rate);
}