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
    *weights = *weights - mul(~(*inputs), *errors)*learning_rate/(double)errors->dims()[0];
    *biases = *biases - (*errors)[0]*learning_rate/(double)errors->dims()[0];
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
    for(size_t i = 1; i < layers.size(); i++)
        layers[layers.size() - i - 1]->backwards(propagator);
}
void Network::predict (Tensor<double>* propagator) {
    for (auto layer : layers) layer->predict(propagator);
}
void Network::fix (double learning_rate) {
    for (auto layer : layers) layer->fix(learning_rate);
}

Model::Model (Layer* network, std::function<double(Tensor<double>*, const Tensor<double>&)> error_func) {
    this->network = network;
    this->error_func = error_func;
}
Model::~Model () {
    delete network;
}
void Model::params (size_t epochs, size_t iterations, size_t batch_size, double learning_rate, double validation_split) {
    this->epochs = epochs;
    this->iterations = iterations;
    this->batch_size = batch_size;
    this->learning_rate = learning_rate;
    this->validation_split = validation_split;
}
void Model::train (const std::unordered_map<std::string, std::vector<double>>& data, const std::vector<std::string>& output_keys) {
    std::unordered_map<std::string, std::vector<double>> data_inputs;
    std::unordered_map<std::string, std::vector<double>> data_outputs;
    std::unordered_map<std::string, std::vector<double>> training_inputs;
    std::unordered_map<std::string, std::vector<double>> training_outputs;
    std::unordered_map<std::string, std::vector<double>> validation_inputs;
    std::unordered_map<std::string, std::vector<double>> validation_outputs;

    for (const auto& pair : data) {
        bool found = false;
        for (const auto& key : output_keys) {
            std::cout << key << " == " << pair.first << '\n';
            if (pair.first.find(key) != std::string::npos) {
                found = true;
                std::cout << "TRUE\n";
                break;
            }
        }
        if (found) data_outputs.insert(pair);
        else data_inputs.insert(pair);
    }
    for (const auto& pair : data_inputs) {
        std::vector<double> training_vector(pair.second.begin(), pair.second.begin() + validation_split*pair.second.size());
        std::vector<double> validation_vector(pair.second.begin() + validation_split*pair.second.size(), pair.second.end());
        training_inputs.insert({pair.first, training_vector});
        validation_inputs.insert({pair.first, validation_vector});
    }
    for (const auto& pair : data_outputs) {
        std::vector<double> training_vector(pair.second.begin(), pair.second.begin() + validation_split*pair.second.size());
        std::vector<double> validation_vector(pair.second.begin() + validation_split*pair.second.size(), pair.second.end());
        training_outputs.insert({pair.first, training_vector});
        validation_outputs.insert({pair.first, validation_vector});
    }
    for (size_t epoch = 0; epoch < epochs; epoch++) {
        size_t output_size;
        size_t input_size;
        for (size_t iteration = 0; iteration < iterations; iteration++) {


            std::vector<double> input_values;
            for (size_t batch = 0; batch < batch_size; batch++) {
                input_size = 0;
                for (const auto& pair : training_inputs) {
                    input_size++;
                    input_values.push_back(pair.second[(batch + iteration*batch_size)%pair.second.size()]);
                }
            }
            std::vector<size_t> input_dims = {batch_size, input_size};
            Tensor<double>* trainer = new Tensor<double>(input_values, input_dims);


            std::vector<double> output_values;
            for (size_t batch = 0; batch < batch_size; batch++) {
                output_size = 0;
                for (const auto& pair : training_outputs) {
                    output_size++;
                    output_values.push_back(pair.second[(batch + iteration*batch_size)%pair.second.size()]);
                }
            }
            std::vector<size_t> output_dims = {batch_size, output_size};
            Tensor<double> expected(output_values, output_dims);


            network->forward(trainer);
            double training_cost = error_func(trainer, expected);
            network->backwards(trainer);
            network->fix(learning_rate);
            delete trainer;
        }


        size_t validation_size = validation_inputs.begin()->second.size();


        std::vector<double> validation_input_values;
        for (size_t i = 0; i < validation_size; i++) {
            for (const auto pair : validation_inputs) {
                validation_input_values.push_back(pair.second[i]);
            }
        }
        std::vector<size_t> validation_input_dims = {validation_size, input_size};
        Tensor<double>* validator = new Tensor<double>(validation_input_values, validation_input_dims);


        std::vector<double> validation_output_values;
        for (size_t i = 0; i < validation_size; i++) {
            for (const auto pair : validation_outputs) {
                validation_output_values.push_back(pair.second[i]);
            }
        }
        std::vector<size_t> validation_output_dims = {validation_size, output_size};
        Tensor<double> validation_expected(validation_output_values, validation_output_dims);


        network->predict(validator);
        double validation_cost = error_func(validator, validation_expected);
        std::cout << "EPOCH " << epoch << ": COST = " << validation_cost << '\n';
        delete validator;

    }
}
void Model::predict (const std::vector<double>& input_values) {

}

void act::sigmoid (Tensor<double>* t) {
    *t = 1/(1 + exp(-(*t)));
}
void act::tanh (Tensor<double>* t) {
    *t = tanh(*t);
}
void act::relu (Tensor<double>* t) {
    *t = max(*t, 0);
}
void act::softmax (Tensor<double>* t) {
    Tensor<double> sums = exp(*t)[1];
    *t = exp(*t)/sums;
}
void act::linear (Tensor<double>* t) {

}
void der::sigmoid (Tensor<double>* t) {
    *t = exp(-(*t))/((1 + exp(-(*t)))*(1 + exp(-(*t))));
}
void der::tanh (Tensor<double>* t) {
    *t = 1 - tanh(*t)*tanh(*t);
}
void der::relu (Tensor<double>* t) {
    *t = (*t > 0);
}
void der::softmax (Tensor<double>* t) {

}
void der::linear (Tensor<double>* t) {

}
double err::regression (Tensor<double>* p, const Tensor<double>& e) {
    *p = *p - e;
    double cost = (abs(*p)[0][1]/(double)(p->dims()[0])).vec()[0];
    return cost;
}
double err::binary (Tensor<double>* p, const Tensor<double>& e) {
    Tensor<double> loss = -(e*log(*p) + ((double)1 - e)*log((double)1 - *p));
    *p = *p - e;
    double cost = (loss[0]/(double)(loss.dims()[0])).vec()[0];
    return cost;
}
double err::categorical (Tensor<double>* p, const Tensor<double>& e) {
    Tensor<double> loss = (-e*log(*p))[1];
    *p = *p - e;
    double cost = (loss[0]/(double)(loss.dims()[0])).vec()[0];
    return cost;
}
double err::multyclass (Tensor<double>* p, const Tensor<double>& e) {
    Tensor<double> loss = (-(e*log(*p) + ((double)1 - e)*log((double)1 - *p)))[1];
    *p = *p - e;
    double cost = (loss[0]/(double)(loss.dims()[0])).vec()[0];
    return cost;
}