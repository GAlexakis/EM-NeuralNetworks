#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <unordered_map>

#include "checker.hpp"

template <class T> class Dictionary {
private:
    std::vector<std::string> keys;
    std::unordered_map<std::string, T> data;
public:
    void push_back(const std::string& key,const  T& value) {
        keys.push_back(key);
        data.insert({key, value});
    }

    void insert(int pos, const std::string& key, const T& value) {
        if (pos >= 0)
            keys.insert(keys.begin() + pos, key);
        else
            keys.insert(keys.end() - pos, key);
        data.insert({key, value});
    }

    T& operator[] (const std::string& key) {
        return data[key];
    }

    T& operator[] (int pos) {
        if (pos >= 0)
            return data[keys[pos]];
        else
            return data[keys[keys.size() - pos]];
    }

    std::string operator() (int pos) const {
        if (pos >= 0)
            return keys[pos];
        else
            return keys[keys.size() - pos];
    }

    size_t size() const {
        return keys.size();
    }

    void print_keys(char seperator = '\t') const {
        for (const auto& key : keys) {
            std::cout << key << seperator;
        }
        std::cout << std::endl;
    }

    std::vector<std::string> get_keys () const {
        return keys;
    }

    std::unordered_map<std::string, T> get_data () const {
        return data;
    }

    void operator= (const Dictionary<T>& src) {
        keys = src.get_keys();
        data = src.get_data();
    }
};

template <class T> Dictionary<std::vector<T>> parse_csv (std::string filepath) {
    Dictionary<std::vector<T>> data;
    std::ifstream file;
    T temp;
    std::string key;
    size_t features = 0;
    std::string count_string;
    file.open(filepath);
    std::getline(file, count_string);
    file.close();
    for (const auto i : count_string) {
        features += i == ',' ? 1 : 0;
        }
    file.open(filepath);
    std::vector<std::string> keys;
    for (size_t i = 0; i < features; i++) {
        std::getline(file, key, ',');
        keys.push_back(key);
        data.push_back(key, {});
    }
    std::getline(file, key, '\n');
    keys.push_back(key);
    data.push_back(key, {});

    bool eof = false;
    while (true) {
        for (size_t i = 0; i < data.size(); i++) {
            if (file >> temp) {
                data[i].push_back(temp);
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

template <class T> Dictionary<std::vector<T>> parse_general_csv (std::string filepath) {
    Dictionary<std::vector<T>> data;
    std::ifstream file;
    std::string key;
    size_t features = 0;
    std::string count_string;
    file.open(filepath);
    std::getline(file, count_string);
    file.close();
    for (const auto i : count_string)
        features += i == ',' ? 1 : 0;

    file.open(filepath);
    for (size_t i = 0; i < features; i++) {
        std::getline(file, key, ',');
        data.push_back(key, {});
    }
    std::getline(file, key, '\n');

    data.push_back(key, {});

    bool eof = false;


    std::vector<std::vector<std::string>> table;
    for (size_t i = 0; i < features + 1; i++)
        table.push_back({});
    while (true) {
        std::vector<std::string> tempValues;
        for (size_t i = 0; i < features; i++) {
            if (std::getline(file, key, ','))
                tempValues.push_back(key);
            else {
                eof = true;
                break;
            }
        }
        if (std::getline(file, key, '\n'))
            tempValues.push_back(key);
        else eof = true;

        if (eof) break;

        size_t index = 0;

        for (size_t j = 0; j < data.size(); j++) {
            if (is_number(tempValues[index]))
                data[j].push_back((T)std::stod(tempValues[index]));
            else {
                bool found = false;
                size_t pos;
                for (size_t i = 0; i < table[index].size(); i++) {
                    if (table[index][i] == tempValues[index]) {
                        found = true;
                        pos = i;
                        break;
                    }
                }
                if (found) {
                    data[j].push_back((T)pos);
                    }
                else {
                    data[j].push_back((T)table[index].size());
                    table[index].push_back(tempValues[index]);
                }
            }
            index++;
        }
    }
    file.close();
    return data;
}

template <class T> void make_one_hot (Dictionary<T>& data, std::string key, std::vector<std::string> one_hot_keys) {
    Dictionary<T> new_data;
    for (int i = 0; i < data.size(); i++){
        if (data(i).find(key) != std::string::npos) {
            for (const auto& one_hot_key : one_hot_keys) {
                new_data.push_back(one_hot_key, {});
            }
            for (int j = 0; j < data[i].size(); j++) {
                for (int k = 0; k < one_hot_keys.size(); k++) {
                    if (data[i][j] == k)
                        new_data[one_hot_keys[k]].push_back(1);
                    else
                        new_data[one_hot_keys[k]].push_back(0);
                }
            }
        }
        else {
            new_data.push_back(data(i), data[i]);
        }
    }
    data =  new_data;
}

template <class T> void print_dictionary (Dictionary<T>& data) {
    data.print_keys();
    for (size_t i = 0; i < data[0].size(); i++) {
        for (size_t j = 0; j < data.size(); j++)
            std::cout << std::setw(data(j).size()) << std::left << data[j][i] << '\t';
        std::cout << '\n';
    }
}