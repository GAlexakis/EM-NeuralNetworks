#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <unordered_map>

template <class T> std::unordered_map<std::string, std::vector<T>> parse_csv (std::string filepath, size_t features) {
    std::unordered_map<std::string, std::vector<float>> data, tempData;
    std::ifstream file;
    T temp;
    std::string key;
    file.open(filepath);
    for (size_t i = 0; i < features - 1; i++) {
        std::getline(file, key, ',');
        tempData.insert({key, {}});
    }
    std::getline(file, key, '\n');
    tempData.insert({key, {}});
    for (const auto &pair : tempData) data.insert(pair);
    bool eof = false;
    while (true) {
        for (auto &pair : data) {
            if (file >> temp) {
                pair.second.push_back(temp);
                if (file.peek() == ',' || file.peek() == '\n' || file.peek() == '\r') file.ignore();
            }
            else {
                eof = true;
                break;
            }
        }
        if (eof) break;
    }
    file.close();
    return data;
}

template <class T> void print_map (const std::unordered_map<std::string,std::vector<T>>& data) {
    size_t size;
    for (const auto &pair : data) {
        std::cout << pair.first << '\t';
        size = pair.second.size();
    }
    std::cout << '\n';
    for (size_t i = 0; i < size; i++) {
        for (const auto &pair : data) std::cout << std::setw(pair.first.size()) << std::left <<pair.second[i] << '\t';
        std::cout << '\n';
    }
}