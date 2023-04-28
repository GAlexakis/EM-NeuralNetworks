#pragma once

#include <functional>

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
    Dense (size_t input_size, size_t output_size, size_t batch_size);
    ~Dense ();
    void forward (Tensor<double>* propagator) override;
    void backwards (Tensor<double>* propagator) override;
    void predict (Tensor<double>* propagator) override;
    void fix (double learning_rate) override;
};

class Activation : public Layer{
private:
    Tensor<double>* inputs;
    std::function<void(Tensor<double>*)> func;
    std::function<void(Tensor<double>*)> der;
public:
    Activation (size_t input_size, size_t batch_size, std::function<void(Tensor<double>*)> func, std::function<void(Tensor<double>*)> der);
    ~Activation ();
    void forward (Tensor<double>* propagator) override;
    void backwards (Tensor<double>* propagator) override;
    void predict (Tensor<double>* propagator) override;
    void fix (double learning_rate) override;
};

class Network : public Layer{
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