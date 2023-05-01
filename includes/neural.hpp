#pragma once

#include <string>
#include <functional>

enum class LAYER {
    DENSE
};
enum class ACTIVATION {
    LINEAR,
    SIGMOID,
    TANH,
    RELU,
    SOFTMAX
};

class Layer {
public:
    virtual ~Layer () {}
    virtual void forward (Tensor<double>* propagator) = 0;
    virtual void backwards (Tensor<double>* propagator) = 0;
    virtual void predict (Tensor<double>* propagator) = 0;
    virtual void fix (double learning_rate) = 0;
};

class Dense : public Layer {
private:
    Tensor<double>* weights;
    Tensor<double>* biases;
    Tensor<double>* inputs;
    Tensor<double>* errors;
public:
    Dense (size_t input_size, size_t output_size);
    ~Dense ();
    void forward (Tensor<double>* propagator) override;
    void backwards (Tensor<double>* propagator) override;
    void predict (Tensor<double>* propagator) override;
    void fix (double learning_rate) override;
};

class Activation : public Layer {
private:
    Tensor<double>* inputs;
    std::function<void(Tensor<double>*)> func;
    std::function<void(Tensor<double>*)> der;
public:
    Activation (size_t input_size, std::function<void(Tensor<double>*)> func, std::function<void(Tensor<double>*)> der);
    ~Activation ();
    void forward (Tensor<double>* propagator) override;
    void backwards (Tensor<double>* propagator) override;
    void predict (Tensor<double>* propagator) override;
    void fix (double learning_rate) override;
};

class Network : public Layer {
private:
    std::vector<Layer*> layers;
public:
    Network ();
    ~Network ();
    void add (Layer* l);
    void forward (Tensor<double>* propagator) override;
    void backwards (Tensor<double>* propagator) override;
    void predict (Tensor<double>* propagator) override;
    void fix (double learning_rate) override;
};

class Model {
private:
    Network* network;
    size_t epochs;
    size_t iterations;
    size_t batch_size;
    size_t last_size;
    double learning_rate;
    double validation_split;
    std::function<double(Tensor<double>*, const Tensor<double>&)> error_func;
public:
    Model (Network* network, std::function<double(Tensor<double>*, const Tensor<double>&)> error_func);
    Model (std::function<double(Tensor<double>*, const Tensor<double>&)> error_func, size_t input_size);
    ~Model ();
    void params (size_t epochs, size_t iterations, size_t batch_size, double learning_rate, double validation_split);
    void train (Dictionary<std::vector<double>>& data, const std::vector<std::string>& output_keys, const std::vector<std::string>& ignore_keys);
    void predict (const std::vector<double>& input_values);
    void add (LAYER l, ACTIVATION a, size_t output_size);
};

namespace act {
    void sigmoid (Tensor<double>* t);
    void tanh (Tensor<double>* t);
    void relu (Tensor<double>* t);
    void softmax (Tensor<double>* t);
    void linear (Tensor<double>* t);
}

namespace der {
    void sigmoid (Tensor<double>* t);
    void tanh (Tensor<double>* t);
    void relu (Tensor<double>* t);
    void softmax (Tensor<double>* t);
    void linear (Tensor<double>* t);
}

namespace err {
    double regression (Tensor<double>* p, const Tensor<double>& e);
    double binary (Tensor<double>* p, const Tensor<double>& e);
    double categorical (Tensor<double>* p, const Tensor<double>& e);
    double multyclass (Tensor<double>* p, const Tensor<double>& e);
}