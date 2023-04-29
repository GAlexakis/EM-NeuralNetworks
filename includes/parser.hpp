#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <unordered_map>

template <class T> std::unordered_map<std::string, std::vector<T>> parse_csv (std::string filepath) {
    std::unordered_map<std::string, std::vector<T>> data, tempData;
    std::ifstream file;
    T temp;
    std::string key;
    size_t features = 0;
    std::string count_string;
    file.open(filepath);
    std::getline(file, count_string);
    file.close();
    for (const auto i : count_string) features += i == ',' ? 1 : 0;
    file.open(filepath);
    for (size_t i = 0; i < features; i++) {
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