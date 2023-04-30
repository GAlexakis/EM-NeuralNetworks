#include<iostream>
#include<vector>
#include<cmath>

#define DEBUG 0
template <typename T> class matrix
{
    T* A;
public:
    size_t rows;
    size_t cols;

    matrix (size_t r = 1, size_t c = 1)
    {
        rows = r;
        cols = c;
        A = new T[rows*cols]{};
        if(DEBUG) std::cout << "creating " << A << '\n';
    }

    matrix (std::vector<T> values, size_t r = 1, size_t c = 1)
    {
        rows = r;
        cols = c;
        A = new T[rows*cols];
        if(DEBUG) std::cout << "creating " << A << '\n';
        for(size_t i = 0; i < rows; i++)
        {
            for(size_t j = 0; j < cols; j++)
            {
                A[i*cols + j] = values[i*cols + j];
            }
        }
    }

    ~matrix ()
    {
        if(DEBUG) std::cout << "deleting " << A << '\n';
        delete[] A;
    }

    T operator() (size_t x , size_t y) const
    {
        return A[x*cols + y];
    }

    void operator() (T value, size_t x, size_t y)
    {
        A[x*cols + y] = value;
    }

    T operator& () const
    {
        T sum = 0;
        for(size_t i = 0; i < rows; i++)
        {
            for(size_t j = 0; j < cols; j++)
            {
                sum += A[i*cols + j];
            }
        }
        return sum;
    }

    matrix<T> operator[] (int s) const
    {
        size_t r, c;
        if(s == 0)
        {
            r = 1;
            c = cols;
        }
        else if (s == 1)
        {
            r = rows;
            c = 1;
        }
        else
        {
            std::cerr << "\x1B[31mERROR: SUM OF ELEMENTS FAILED\033[0m\n";
            exit(1);
        }

        matrix<T> result(r, c);
        if(s == 0)
        {
            for(size_t i = 0; i < cols; i++)
            {
                result(0, 0, i);
                for(size_t j = 0; j < rows; j++)
                {
                    result(result(0, i) + A[i + j*cols], 0, i);
                }
            }
        }
        else if (s == 1)
        {
            for(size_t i = 0; i < rows; i++)
            {
                result(0, i, 0);
                for(size_t j = 0; j < cols; j++)
                {
                    result(result(i, 0) + A[i*cols + j], i, 0);
                }
            }
        }
        return result;
    }

    void operator<< (char seperator) const
    {
        for(size_t i = 0; i < rows; i++)
        {
            for(size_t j = 0; j < cols; j++)
            {
                std::cout << A[i*cols + j] << seperator;
            }
            std::cout << '\n';
        }
    }

    void operator>> (const matrix<T>& M)
    {
        if(DEBUG) std::cout << "reallocating " << A << '\n';
        delete[] A;
        rows = M.rows;
        cols = M.cols;
        A = new T[rows*cols];
        for(size_t i = 0; i < rows; i++)
        {
            for(size_t j = 0; j < cols; j++)
            {
                A[i*cols + j] = M(i, j);
            }
        }
    }

    void operator<= (const matrix<T>& M)
    {
        if(rows != M.rows || cols != M.cols)
        {
            std::cerr << "\x1B[31mERROR: COPY FAILED\033[0m\n";
            exit(1);
        }
        for(size_t i = 0; i < rows; i++)
        {
            for(size_t j = 0; j < cols; j++)
            {
                A[i*cols + j] = M(i, j);
            }
        }
    }

    matrix<T> operator* (T a) const
    {
        matrix<T> result(rows, cols);
        for(size_t i = 0; i < rows; i++)
        {
            for(size_t j = 0; j < cols; j++)
            {
                result(a*A[i*cols + j], i, j);
            }
        }


        return result;
    }

    matrix<T> operator* (const matrix<T>& M) const
    {
        matrix<T> result(rows, cols);
        if(M.rows == rows && M.cols == cols)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j]*M(i, j), i, j);
                }
            }
        }
        else if(M.rows == rows && M.cols == 1)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j]*M(i, 0), i, j);
                }
            }
        }
        else if(M.rows == 1 && M.cols == cols)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j]*M(0, j), i, j);
                }
            }
        }
        else
        {
            std::cerr << "\x1B[31mERROR: ELEMENT WISE MULTIPLICATION FAILED\033[0m\n";
            exit(1);
        }
        return result;
    }

    matrix<T> operator/ (T a) const
    {
        matrix<T> result(rows, cols);
        for(size_t i = 0; i < rows; i++)
        {
            for(size_t j = 0; j < cols; j++)
            {
                result(A[i*cols + j]/a, i, j);
            }
        }
        return result;
    }

    matrix<T> operator/ (const matrix<T>& M) const
    {
        matrix<T> result(rows, cols);
        if(M.rows == rows && M.cols == cols)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j]/M(i, j), i, j);
                }
            }
        }
        else if(M.rows == rows && M.cols == 1)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j]/M(i, 0), i, j);
                }
            }
        }
        else if(M.rows == 1 && M.cols == cols)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j]/M(0, j), i, j);
                }
            }
        }
        else
        {
            std::cerr << "\x1B[31mERROR: ELEMENT WISE DIVISION FAILED\033[0m\n";
            exit(1);
        }
        return result;
    }

    matrix<T> operator+ (T a) const
    {
        matrix<T> result(rows, cols);
        for(size_t i = 0; i < result.rows; i++)
        {
            for(size_t j = 0; j < result.cols; j++)
            {
                result(A[i*cols + j] + a, i, j);
            }
        }
        return result;
    }

    matrix<T> operator+ (const matrix<T>& M) const
    {
        matrix<T> result(rows, cols);
        if(M.rows == rows && M.cols == cols)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j] + M(i, j), i, j);
                }
            }
        }
        else if(M.rows == rows && M.cols == 1)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j] + M(i, 0), i, j);
                }
            }
        }
        else if(M.rows == 1 && M.cols == cols)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j] + M(0, j), i, j);
                }
            }
        }
        else
        {
            std::cerr << "\x1B[31mERROR: MATRIX ADDITION FAILED\033[0m\n";
            exit(1);
        }
        return result;
    }

    matrix<T> operator- (T a) const
    {
        matrix<T> result(rows, cols);
        for(size_t i = 0; i < result.rows; i++)
        {
            for(size_t j = 0; j < result.cols; j++)
            {
                result(A[i*cols + j] - a, i, j);
            }
        }
        return result;
    }

    matrix<T> operator- (const matrix<T>& M) const
    {
        matrix<T> result(rows, cols);
        if(M.rows == rows && M.cols == cols)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j] - M(i, j), i, j);
                }
            }
        }
        else if(M.rows == rows && M.cols == 1)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j] - M(i, 0), i, j);
                }
            }
        }
        else if(M.rows == 1 && M.cols == cols)
        {
            for(size_t i = 0; i < rows; i++)
            {
                for(size_t j = 0; j < cols; j++)
                {
                    result(A[i*cols + j] - M(0, j), i, j);
                }
            }
        }
        else
        {
            std::cerr << "\x1B[31mERROR: MATRIX SUBSTRACTION FAILED\033[0m\n";
            exit(1);
        }
        return result;
    }

    matrix<T> operator- () const
    {
        matrix<T> result(rows, cols);
        for(size_t i = 0; i < result.rows; i++)
        {
            for(size_t j = 0; j < result.cols; j++)
            {
                result(-A[i*cols + j], i, j);
            }
        }
        return result;
    }

    matrix<T> operator~ ()
    {
        matrix<T> result(cols, rows);
        for(size_t i = 0; i < result.rows; i++)
        {
            for(size_t j = 0; j < result.cols; j++)
            {
                result(A[j*cols + i], i, j);
            }
        }
        return result;
    }

    matrix<T> operator! () const
    {
        matrix<T> result(rows, cols);
        for(size_t i = 0; i < result.rows; i++)
        {
            for(size_t j = 0; j < result.cols; j++)
            {
                result(1/A[i*cols + j], i, j);
            }
        }
        return result;
    }

    matrix<T> operator% (const matrix<T>& M) const
    {
        if(cols != M.rows)
        {
            std::cerr << "\x1B[31mERROR: MATRIX MULTIPLICATION FAILED\033[0m\n";
            exit(1);
        }
        matrix<T> result(rows, M.cols);
        for(size_t i = 0; i < result.rows; i++)
        {
            for(size_t j = 0; j < result.cols; j++)
            {
                result(0, i, j);
                for(size_t k = 0; k < cols; k++)
                {
                    result(result(i, j) + A[i*cols + k]*M(k,j), i, j);
                }
            }
        }
        return result;
    }
};

typedef enum LayerType
{
    DENSE,
    ACTIVATION
}
layerType_t;

typedef enum Activation
{
    SIGMOID,
    RELU,
    TANH,
    LINEAR,
    SOFTMAX,
}
activation_t;

matrix<double> exp (const matrix<double>& M)
{
    matrix<double> result(M.rows, M.cols);
    for(size_t i = 0; i < M.rows; i++)
    {
        for(size_t j = 0; j < M.cols; j++)
        {
            result(exp(M(i, j)), i, j);
        }
    }
    return result;
}

matrix<double> tanh (const matrix<double>& M)
{
    matrix<double> result(M.rows, M.cols);
    for(size_t i = 0; i < M.rows; i++)
    {
        for(size_t j = 0; j < M.cols; j++)
        {
            result(tanh(M(i, j)), i, j);
        }
    }
    return result;
}

matrix<double> max (const matrix<double>& M, double a)
{
    matrix<double> result(M.rows, M.cols);
    for(size_t i = 0; i < M.rows; i++)
    {
        for(size_t j = 0; j < M.cols; j++)
        {
            result(fmax(M(i, j), a), i, j);
        }
    }
    return result;
}

matrix<double> abs (const matrix<double>& M)
{
    matrix<double> result(M.rows, M.cols);
    for(size_t i = 0; i < M.rows; i++)
    {
        for(size_t j = 0; j < M.cols; j++)
        {
            result(abs(M(i, j)), i, j);
        }
    }
    return result;
}

matrix<double> sigmoid (const matrix<double>& M)
{
    return !(exp(-M) + 1);
}

matrix<double> relu (const matrix<double>& M)
{
    return max(M, 0);
}

matrix<double> softmax (const matrix<double>& M)
{
    matrix<double> result = exp(M);
    double sum = &result;
    std::cout << sum << '\n';
    result <= result/sum;
    return result;
}

matrix<double> sigmoid_der (const matrix<double>& M)
{
    return sigmoid(M)*(-sigmoid(M) + 1);
}

matrix<double> tanh_der (const matrix<double>& M)
{
    return -(tanh(M)*tanh(M)) + 1;
}

matrix<double> relu_der (const matrix<double>& M)
{
    matrix<double> result(M.rows, M.cols);
    for(size_t i = 0; i < M.rows; i++)
    {
        for(size_t j = 0; j < M.cols; j++)
        {
            result((M(i, j) > 0)? 1 : 0, i, j);
        }
    }
    return result;
}

matrix<double> softmax_der (const matrix<double>& M)
{
    return softmax(M)*(-softmax(M) + 1);
}

matrix<double>* regression_cost (const matrix<double>& predictions, const matrix<double>& labels, size_t batch_size)
{
    matrix<double>* result = new matrix<double> (1,1);
    *result >> abs(predictions - labels)[0]/(double)batch_size;
    // *result << ',';
    return result;
}

void regression_error (matrix<double>* propagator, const matrix<double>& labels)
{
    (*propagator) <= ((*propagator) - labels);
    // *propagator << ',';
}

typedef class Layer
{
    activation_t activation;
    layerType_t layerType;
    size_t input_size;
    size_t output_size;
    size_t batch_size;

    matrix<double> *weights;
    matrix<double> *biases;
    matrix<double> *inputs;
    matrix<double> *outputs;
    matrix<double> *errors;

public:
    Layer (size_t in, size_t out, layerType_t lrt = DENSE, activation_t act = LINEAR, size_t bs = 1)
    {
        input_size = in;
        output_size = out;
        batch_size = bs;
        layerType = lrt;
        activation = act;

        weights = new matrix<double>(input_size, output_size);
        biases = new matrix<double>(1, output_size);
        inputs = new matrix<double>(batch_size, input_size);
        outputs = new matrix<double>(batch_size, output_size);
        errors = new matrix<double>(batch_size, output_size);
        for( size_t i = 0; i < input_size; i++)
        {
            for( size_t j = 0; j < output_size; j++)
            {
                (*weights)(rand()/(double)(RAND_MAX) - 0.5, i, j);
            }
        }
        for( size_t i = 0; i < output_size; i++)
        {
            (*biases)(0, 0, i);
        }
    }

    ~Layer ()
    {
        delete weights;
        delete biases;
        delete inputs;
        delete outputs;
        delete errors;
    }

    void forward (matrix<double>* input)
    {
        // std::cout << "FORWARD PROPAGATION\n";
        (*inputs) <= (*input);
        // (*input) << ',';
        // std::cout << "*\n";
        // (*weights) << ',';
        // std::cout << "+\n";
        // (*biases) << ',';
        // std::cout << "=\n";
        (*input) >> (*input)%(*weights) + (*biases);
        // (*input) << ',';
        (*outputs) <= (*input);
        switch(activation)
        {
        case SIGMOID:
            // std::cout << "~(sigmoid)~\n";
            (*input) <= sigmoid(*input);
            // (*input) << ',';
            break;

        case RELU:
            // std::cout << "~(relU)~\n";
            (*input) <= relu(*input);
            // (*input) << ',';
            break;

        case TANH:
            // std::cout << "~(tanh)~\n";
            (*input) <= tanh(*input);
            // (*input) << ',';
            break;

        case SOFTMAX:
            // std::cout << "~(softmax)~\n";
            (*input) <= softmax(*input);
            // (*input) << ',';
            break;

        default:
            break;
        }
        // std::cout << "----------------------------------------\n";
    }

    void backwards (matrix<double>* delta)
    {
        // std::cout << "BACKWARDS PROPAGATION\n";
        // *delta << ',';

        switch(activation)
        {
        case SIGMOID:
            // std::cout << "* sigmoid der\n";
            *delta <= (*delta)*sigmoid_der(*outputs);
            // *delta << ',';
            break;

        case RELU:
            // std::cout << "* relU der\n";
            *delta <= (*delta)*relu_der(*outputs);
            // *delta << ',';
            break;

        case TANH:
            // std::cout << "* tanh der\n";
            *delta <= (*delta)*tanh_der(*outputs);
            // *delta << ',';
            break;

        case SOFTMAX:
            // std::cout << "* softmax der\n";
            *delta <= (*delta)*softmax_der(*outputs);
            // *delta << ',';
            break;

        default:
            break;
        }

        (*errors) <= (*delta);
        // std::cout << "* Transpose\n";
        // *weights << ',';

        // std::cout << "=\n";

        *delta >> (*delta)%(~(*weights));

        // *delta << ',';

        // std::cout << "----------------------------------------\n";
    }

    void update (double learning_rate = 0.01)
    {
        // std::cout << "UPDATING PARAMETERS\n";
        // (*weights) << ',';
        // std::cout << "||||||\n";
        (*weights) <= (*weights) - ((~(*inputs))%(*errors))*learning_rate/(double)batch_size;
        // (*weights) << ',';
        // std::cout << "\n";
        // (*biases) << ',';
        // std::cout << "-\n";
        // (*errors) << ',';
        // std::cout << "*" << learning_rate << '\n';
        // std::cout << "=\n";
        (*biases) <= (*biases) - (*errors)[0]*learning_rate/(double)batch_size;
        // (*biases) << ',';
        // std::cout << "----------------------------------------\n";
    }
}
layer_t;

typedef class Network{
    size_t input_size;
    size_t output_size;
    size_t split;
    size_t batch_size;
    std::vector<double> samples;
    std::vector<double> labels;
    std::vector<layer_t*> layers;

public:
    Network (std::vector<double> s, std::vector<double> l, size_t ins, size_t outs, size_t sp)
    {
        input_size = ins;
        output_size = outs;
        split = sp;
        samples = s;
        labels = l;
        layers = {};
        batch_size = (size_t)(samples.size()/(split*input_size));
    }

    ~Network ()
    {
        for(size_t i = 0; i < layers.size(); i++)
        {
            delete layers[i];
        }
    }

    void create (size_t in, size_t out, layerType_t lrt = DENSE, activation_t act = LINEAR)
    {
        layer_t* lptr = new layer_t(in, out, lrt, act, batch_size);
        layers.push_back(lptr);
    }

    void train(size_t epochs, double learning_rate)
    {
        matrix<double>* loss;
        matrix<double>* propagator;
        matrix<double>* evaluator;
        for(size_t i = 0; i < epochs; i++)
        {
            double cost = 0;
            std::cout << "epoch: " << i << '\n';
            for(size_t j = 0; j < split; j++)
            {

                std::vector<double> sample_batch(   samples.begin() + input_size*batch_size*j*samples.size(),
                                                    samples.begin() + input_size*batch_size*(j + 1)*samples.size());
                propagator = new matrix<double>(sample_batch, batch_size, input_size);

                for(size_t k = 0; k < layers.size(); k++)
                {
                    layers[k]->forward(propagator);
                }

                std::vector<double> label_batch(    labels.begin() + output_size*batch_size*j*labels.size(),
                                                    labels.begin() + output_size*batch_size*(j + 1)*labels.size());
                evaluator = new matrix<double>(label_batch, batch_size, output_size);

                loss = regression_cost(*propagator, *evaluator, batch_size);
                regression_error(propagator, *evaluator);

                for(size_t k = 0; k < layers.size(); k++)
                {
                    layers[layers.size() - k - 1]->backwards(propagator);
                }
                delete propagator;
                for(size_t k = 0; k < layers.size(); k++)
                {
                    layers[k]->update(learning_rate);
                }
                cost += &(*loss);
                delete loss;
            }
            std::cout << "cost: " << cost << '\n';
        }
    }
}
network_t;

int main (int argc, char** argv)
{
    network_t net(
        {
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10
        },
        {
            5*1,
            5*2,
            5*3,
            5*4,
            5*5,
            5*6,
            5*7,
            5*8,
            5*9,
            5*10
        },1, 1, 1
    );
    net.create(1, 1, DENSE, LINEAR);
    // net.create(1, 1, DENSE, LINEAR);

    net.train(100, 0.01);
    return 0;
}