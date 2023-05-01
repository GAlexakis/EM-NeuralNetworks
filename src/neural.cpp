#include "parser.hpp"
#include "tensor.hpp"
#include "progress.hpp"
#include "neural.hpp"

Dense::Dense (size_t input_size, size_t output_size) {
    weights = new Tensor<double>(input_size, output_size);
    biases  = new Tensor<double>((size_t)1, output_size);
    inputs  = new Tensor<double>;
    errors  = new Tensor<double>;

    for(size_t i = 0; i < input_size; i++) {
            for( size_t j = 0; j < output_size; j++) {
                (*weights)(i, j) = rand()/(double)(RAND_MAX) - 0.5;
            }
        }
        for(size_t i = 0; i < output_size; i++) {
            (*biases)((size_t)0, i) = 0;
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

Activation::Activation (size_t input_size, std::function<void(Tensor<double>*)> func, std::function<void(Tensor<double>*)> der) {
    inputs  = new Tensor<double>;
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

Model::Model (Network* network, std::function<double(Tensor<double>*, const Tensor<double>&)> error_func) {
    this->network = network;
    this->error_func = error_func;
    this->last_size = 0;
}
Model::Model (std::function<double(Tensor<double>*, const Tensor<double>&)> error_func, size_t input_size) {
    this->network = new Network();
    this->error_func = error_func;
    this->last_size = input_size;
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
void Model::train (Dictionary<std::vector<double>>& data, const std::vector<std::string>& output_keys, const std::vector<std::string>& ignore_keys) {
    Dictionary<std::vector<double>> data_inputs;
    Dictionary<std::vector<double>> data_outputs;
    Dictionary<std::vector<double>> training_inputs;
    Dictionary<std::vector<double>> training_outputs;
    Dictionary<std::vector<double>> validation_inputs;
    Dictionary<std::vector<double>> validation_outputs;


    for (size_t i = 0; i < data.size(); i++) {
        bool ignorefound = false;
        for (const auto& key : ignore_keys) {
            if (data(i).find(key) != std::string::npos) {
                ignorefound = true;
                break;
            }
        }
        if (ignorefound) continue;
        bool found = false;
        for (const auto& key : output_keys) {
            if (data(i).find(key) != std::string::npos) {
                found = true;
                break;
            }
        }
        if (found) data_outputs.push_back(data(i), data[i]);
        else data_inputs.push_back(data(i), data[i]);
    }
    for (size_t i = 0; i < data_inputs.size(); i++) {
        std::vector<double> training_vector(data_inputs[i].begin(), data_inputs[i].begin() + validation_split*data_inputs[i].size());
        std::vector<double> validation_vector(data_inputs[i].begin() + validation_split*data_inputs[i].size(), data_inputs[i].end());
        training_inputs.push_back(data_inputs(i), training_vector);
        validation_inputs.push_back(data_inputs(i), validation_vector);
    }
    for (size_t i = 0; i < data_outputs.size(); i++) {
        std::vector<double> training_vector(data_outputs[i].begin(), data_outputs[i].begin() + validation_split*data_outputs[i].size());
        std::vector<double> validation_vector(data_outputs[i].begin() + validation_split*data_outputs[i].size(), data_outputs[i].end());
        training_outputs.push_back(data_outputs(i), training_vector);
        validation_outputs.push_back(data_outputs(i), validation_vector);
    }
    for (size_t epoch = 0; epoch < epochs; epoch++) {
        Progress bar({(int)iterations});
        bar.show({""});
        size_t output_size;
        size_t input_size;
        double average_cost;
        double sum_cost = 0;
        for (size_t iteration = 0; iteration < iterations; iteration++) {


            std::vector<double> input_values;
            for (size_t batch = 0; batch < batch_size; batch++) {
                input_size = 0;
                for (size_t i = 0; i < training_inputs.size(); i++) {
                    input_size++;
                    input_values.push_back(training_inputs[i][(batch + iteration*batch_size)%training_inputs[i].size()]);
                }
            }
            std::vector<size_t> input_dims = {batch_size, input_size};
            Tensor<double>* trainer = new Tensor<double>(input_values, input_dims);


            std::vector<double> output_values;
            for (size_t batch = 0; batch < batch_size; batch++) {
                output_size = 0;
                for (size_t i = 0; i < training_outputs.size(); i++) {
                    output_size++;
                    output_values.push_back(training_outputs[i][(batch + iteration*batch_size)%training_outputs[i].size()]);
                }
            }
            std::vector<size_t> output_dims = {batch_size, output_size};
            Tensor<double> expected(output_values, output_dims);


            network->forward(trainer);
            double training_cost = error_func(trainer, expected);
            network->backwards(trainer);
            network->fix(learning_rate);
            delete trainer;

            sum_cost += training_cost;
            average_cost = sum_cost/(iteration + 1);

            bar.update(0);
            bar.show({" EPOCH: "+ std::to_string(epoch) + " TRAINING COST: " + std::to_string(average_cost) + "   "});

        }


        size_t validation_size = validation_inputs[0].size();


        std::vector<double> validation_input_values;
        for (size_t i = 0; i < validation_size; i++) {
            for (size_t j = 0; j < validation_inputs.size(); j++) {
                validation_input_values.push_back(validation_inputs[j][i]);
            }
        }
        std::vector<size_t> validation_input_dims = {validation_size, input_size};
        Tensor<double>* validator = new Tensor<double>(validation_input_values, validation_input_dims);


        std::vector<double> validation_output_values;
        for (size_t i = 0; i < validation_size; i++) {
            for (size_t j = 0; j < validation_outputs.size(); j++) {
                validation_output_values.push_back(validation_outputs[j][i]);
            }
        }
        std::vector<size_t> validation_output_dims = {validation_size, output_size};
        Tensor<double> validation_expected(validation_output_values, validation_output_dims);


        network->predict(validator);
        double validation_cost = error_func(validator, validation_expected);
        delete validator;

        bar.show({" EPOCH: "+ std::to_string(epoch) + " TRAINING COST: " + std::to_string(average_cost) + " VALIDATION COST: " + std::to_string(validation_cost) + "   "});

    }
}
void Model::predict (const std::vector<double>& input_values) {
    Tensor<double>* predictor = new Tensor<double>(input_values, {1,input_values.size()});
    network->predict(predictor);
    predictor->print();
    delete predictor;
}
void Model::add (LAYER l, ACTIVATION a, size_t output_size) {
    switch(l) {
    case LAYER::DENSE:
        network->add(new Dense(last_size, output_size));
        last_size = output_size;
        break;
    default:
        break;
    }

    switch(a) {
        case ACTIVATION::LINEAR:
            network->add(new Activation(last_size, act::linear, der::linear));
            break;
        case ACTIVATION::SIGMOID:
            network->add(new Activation(last_size, act::sigmoid, der::sigmoid));
            break;
        case ACTIVATION::TANH:
            network->add(new Activation(last_size, act::tanh, der::tanh));
            break;
        case ACTIVATION::RELU:
            network->add(new Activation(last_size, act::relu, der::relu));
            break;
        case ACTIVATION::SOFTMAX:
            network->add(new Activation(last_size, act::softmax, der::softmax));
            break;
        default:
        break;
    }
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