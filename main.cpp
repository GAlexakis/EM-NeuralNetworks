#include "parser.hpp"
#include "tensor.hpp"
#include "neural.hpp"
#include "progress.hpp"

int main (int argc, char** argv) {

    signal(SIGINT, restore_terminal);

    int example;
    std::cout << "CHOOSE AN EXAMPLE:\n\nDIABETES \t(BINARY) [1]\nHEART \t\t(BINARY) [2]\nF = 5*X + 10*Y \t(REGRESSION) [3]\nY = 3*X \t(REGRESSION) [4]\nSTELLAR \t(CATEGORICAL) [5]\n IRIS \t\t(CATEGORICAL) [6]\n";
    std::cin >> example;
    switch (example) {
    case 1: {
        Dictionary<std::vector<double>> data;
        data = parse_csv<double>("datasets/diabetes.csv");

        std::cout << "DATASET:";
        std::cin.get();

        print_dictionary(data);

        std::cout << "TRAIN NETWORK:";
        std::cin.get();

        Model mod(err::binary, 8);
        mod.add(LAYER::DENSE, ACTIVATION::RELU, 4);
        mod.add(LAYER::DENSE, ACTIVATION::RELU, 4);
        mod.add(LAYER::DENSE, ACTIVATION::SIGMOID, 1);
        mod.params(400, 6 ,100, 0.01, 0.8);
        mod.train(data, {"Outcome"}, {});
    }
        break;
    case 2: {
        Dictionary<std::vector<double>> data;
        data = parse_csv<double>("datasets/heart.csv");

        std::cout << "DATASET:";
        std::cin.get();

        print_dictionary(data);

        std::cout << "TRAIN NETWORK:";
        std::cin.get();

        Model mod(err::binary, 13);
        mod.add(LAYER::DENSE, ACTIVATION::TANH, 8);
        mod.add(LAYER::DENSE, ACTIVATION::TANH, 8);
        mod.add(LAYER::DENSE, ACTIVATION::SIGMOID, 1);
        mod.params(400, 4,60, 0.001, 0.8);
        mod.train(data, {"output"}, {});
    }
        break;
    case 3: {
        Dictionary<std::vector<double>> data;
        data.push_back("val1", {});
        data.push_back("val2", {});
        data.push_back("res", {});
        for (int i = 0; i < 10000; i++) {
           double val1 = (i/100)%100;
           double val2 = i%100;
           double res = 5*val1 + 10*val2;

           data["val1"].push_back(val1);
           data["val2"].push_back(val2);
           data["res"].push_back(res);
        }

        std::cout << "DATASET:";
        std::cin.get();

        print_dictionary(data);

        std::cout << "TRAIN NETWORK:";
        std::cin.get();

        Model mod(err::regression, 2);
        mod.add(LAYER::DENSE, ACTIVATION::LINEAR, 1);
        mod.params(100, 100, 100 , 0.0001, 0.8);
        mod.train(data, {"res"}, {});

        while (true) {
            double val1, val2;
            std::cout << "X = ";
            std::cin >> val1;
            std::cout << "Y = ";
            std::cin >> val2;
            mod.predict({val1, val2});
            std::cout << "REAL VALUE: " << 5*val1 + 10*val2 << '\n';
            if (val1 == 0 && val2 == 0) break;
        }
    }
        break;
    case 4: {
        Dictionary<std::vector<double>> data;
        data.push_back("val1", {});
        data.push_back("res", {});
        for (int i = 0; i < 100; i++) {
           double val1 = i;
           double res = 3*val1;

           data["val1"].push_back(val1);
           data["res"].push_back(res);
        }

        std::cout << "DATASET:";
        std::cin.get();

        print_dictionary(data);

        std::cout << "TRAIN NETWORK:";
        std::cin.get();

        Model mod(err::regression, 1);
        mod.add(LAYER::DENSE, ACTIVATION::LINEAR, 1);
        mod.params(4000, 1, 100 , 0.001, 0.8);
        mod.train(data, {"res"}, {});

        while (true) {
            double val1;
            std::cout << "X = ";
            std::cin >> val1;
            mod.predict({val1});
            std::cout << "REAL VALUE: " << 3*val1 << '\n';
            if (val1 == 0) break;
        }
    }
        break;
    case 5: {
        Dictionary<std::vector<double>> data;
        data = parse_general_csv<double>("datasets/stellar.csv");

        std::cout << "DATASET:";
        std::cin.get();

        print_dictionary(data);

        std::cout << "ONE HOT ENCODEING FOR OUTPUT:";
        std::cin.get();

        make_one_hot(data, "class", {"galaxy", "QSQ", "star"});
        print_dictionary(data);

        std::cout << "TRAIN NETWORK:";
        std::cin.get();

        Model mod(err::categorical, 11);
        mod.add(LAYER::DENSE, ACTIVATION::TANH, 8);
        mod.add(LAYER::DENSE, ACTIVATION::TANH, 8);
        mod.add(LAYER::DENSE, ACTIVATION::SOFTMAX, 3);
        mod.params(40,100,80,0.01,0.8);
        mod.train(data, {"galaxy", "QSQ", "star"}, {"obj_ID", "run_ID", "rerun_ID", "field_ID", "spec_obj_ID", "fiber_ID"});
    }
        break;
    case 6: {
        Dictionary<std::vector<double>> data;
        data = parse_general_csv<double>("datasets/iris.csv");

        std::cout << "DATASET:";
        std::cin.get();

        print_dictionary(data);

        std::cout << "ONE HOT ENCODEING FOR OUTPUT:";
        std::cin.get();

        make_one_hot(data, "species", {"setosa", "versicolor", "virginica"});
        print_dictionary(data);

        std::cout << "TRAIN NETWORK:";
        std::cin.get();

        Model mod(err::categorical, 4);
        mod.add(LAYER::DENSE, ACTIVATION::TANH, 4);
        mod.add(LAYER::DENSE, ACTIVATION::SOFTMAX, 3);
        mod.params(4000,1,120,0.01,0.8);
        mod.train(data, {"setosa", "versicolor", "virginica"}, {});
        while (true) {
            double val1, val2, val3, val4;
            std::cout << "Sepal length = ";
            std::cin >> val1;
            if (std::cin.peek() == ',') std::cin.ignore();
            std::cout << "Sepal width = ";
            std::cin >> val2;
            if (std::cin.peek() == ',') std::cin.ignore();
            std::cout << "Petal length = ";
            std::cin >> val3;
            if (std::cin.peek() == ',') std::cin.ignore();
            std::cout << "Petal width = ";
            std::cin >> val4;
            mod.predict({val1, val2, val3, val4});
            if (val1 == 0 && val2 == 0 && val3 == 0 && val4 == 0) break;
        }
    }
        break;
    default:
        break;
    }

    std::cout << "EXITING PROGRAM\n";
    return 0;
}
