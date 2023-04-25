#include "parser.hpp"

int main (int argc, char** argv) {
    std::unordered_map<std::string, std::vector<float>> data = parse_csv<float>("diabetes.csv", 9);
    print_map(data);
    return 0;
}