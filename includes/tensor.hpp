#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <omp.h>

template<class T> class Tensor {
private:
    size_t size;
    std::vector<size_t> dimensions;
    T* ptr;

public:
    Tensor (std::vector<size_t> dims) {
        size = 1;
        dimensions = dims;
        for (const size_t& dimension : dimensions) {
            size *= dimension;
        }
        ptr = new T[size]{};
    }

    template <class... Args> Tensor (Args ...args) {
        size = 1;
        dimensions = {args...};
        for (const size_t& dimension : dimensions) {
            size *= dimension;
        }
        ptr = new T[size]{};
    }

    Tensor (const std::vector<T>& values, const std::vector<size_t>& dims) {
        size = 1;
        dimensions = dims;
        for (const size_t& dimension : dimensions) {
            size *= dimension;
        }
        if (values.size() != size) {
            std::cerr << "\x1B[31mERROR: TENSOR SIZE AND ELEMENTS DON'T MATCH\033[0m\n";
            exit(1);
        }
        ptr = new T[size]{};
        for(size_t i = 0; i < size; i++) {
            std::vector<size_t> indices(dimensions.size());
            size_t divider = size/dimensions[0];
            size_t index = 0;
            for (size_t j = 0; j < indices.size(); j++) {
                indices[j] = (i/divider)%dimensions[j];
                index += divider*indices[j];
                divider /= j == indices.size() - 1 ? 1 : dimensions[j + 1];
            }
            ptr[index] = values[index];
        }
    }

    Tensor (const Tensor<T>& src) {
        size = src.length();
        dimensions = src.dims();
        ptr = new T[size];
        for(size_t i = 0; i < size; i++)
        {
            std::vector<size_t> indices(dimensions.size());
            size_t divider = size/dimensions[0];
            size_t index = 0;
            for (size_t j = 0; j < indices.size(); j++) {
                indices[j] = (i/divider)%dimensions[j];
                index += divider*indices[j];
                divider /= j == indices.size() - 1 ? 1 : dimensions[j + 1];
            }
            ptr[index] = src(indices);
        }
    }

    Tensor (Tensor<T>&& src) {
        size = src.length();
        dimensions = src.dims();
        ptr = &src;
        src.del();
    }

    void print(char seperator = ',') const {
        std::cout << '\n';
        std::vector<size_t> temp(dimensions.size());
        std::vector<size_t> indices(dimensions.size());
        size_t dimsize = dimensions.size();
        size_t changed = dimsize - 1;
        for(size_t i = 0; i < size; i++) {
            size_t divider = size/dimensions[0];
            size_t index = 0;
            bool found = false;
            for (size_t j = 0; j < indices.size(); j++) {
                temp[j] = (i/divider)%dimensions[j];
                if (i != 0 && !found && indices[j] != temp[j]) {
                    changed = j;
                    found = true;
                }
                indices[j] = temp[j];
                index += divider*indices[j];
                divider /= j == indices.size() - 1 ? 1 : dimensions[j + 1];
            }
            if (i != 0) {
                if (changed == dimsize - 1) std::cout << seperator;
                else if (changed == dimsize - 2) std::cout << '\n';
                else if (changed == dimsize - 3) std::cout << "\n\n";
                else if (changed == dimsize - 4) std::cout << "\n---------------------------------------------------------------------------------------------------- DIMENSION " << changed + 1 <<"\n";
                else if (changed == dimsize - 5) std::cout << "\n==================================================================================================== DIMENSION " << changed + 1 <<"\n";
                else if (changed == dimsize - 6) std::cout << "\n----------------------------------------------------------------------------------------------------\n|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| DIMENSION " << changed + 1 <<"\n----------------------------------------------------------------------------------------------------\n";
                else if (changed == dimsize - 7) std::cout << "\n====================================================================================================\n|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| DIMENSION " << changed + 1 <<"\n====================================================================================================\n";
                else std::cout << "\n====================================================================================================\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n---------------------------------------------------------------------------------------------------- DIMENSION " << changed + 1 <<"\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n====================================================================================================\n";
            }
            std::cout << ptr[index];
        }
        std::cout << '\n';
    }

    std::vector<T> vec () {
        std::vector<T> result(size);
        for(size_t i = 0; i < size; i++) {
            result[i] = ptr[i];
        }
        return result;
    }

    ~Tensor () {
        delete[] ptr;
    }

    void del () {
        size = 0;
        dimensions = {};
        ptr = nullptr;
    }

    size_t length() const {
        return size;
    }

    std::vector<size_t> dims () const {
        return dimensions;
    }

    T* operator& () const {
        return ptr;
    }

    template <class... Args> T& operator() (Args ...args) const {
        std::vector<int> indices = {args...};
        if (indices.size() != dimensions.size()) {
            std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
            exit(1);
        }
        size_t divider = size/dimensions[0];
        size_t index = 0;
        for (size_t i = 0; i < indices.size(); i++) {
            index += divider*indices[i];
            divider /= i == indices.size() - 1 ? 1 : dimensions[i + 1];
            if (index >= size) {
                std::cerr << "\x1B[31mERROR: OUT OF BOUNDS\033[0m\n";
                exit(1);
            }
        }
        return ptr[index];
    }

    T& operator() (std::vector<size_t> indices) const {
        if (indices.size() != dimensions.size()) {
            std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
            exit(1);
        }
        size_t divider = size/dimensions[0];
        size_t index = 0;
        for (size_t i = 0; i < indices.size(); i++) {
            index += divider*indices[i];
            divider /= i == indices.size() - 1 ? 1 : dimensions[i + 1];
            if (index >= size) {
                std::cerr << "\x1B[31mERROR: OUT OF BOUNDS\033[0m\n";
                exit(1);
            }
        }
        return ptr[index];
    }

    Tensor<T> operator[] (size_t dimension) const {
        if (dimension >= dimensions.size() || dimension < 0) {
            std::cerr << "\x1B[31mERROR: OUT OF DIMENSION BOUNDS\033[0m\n";
            exit(1);
        }
        std::vector<size_t> new_dims = dimensions;
        new_dims[dimension] = 1;
        Tensor<T> result = Tensor<T>(new_dims);
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(result.dims().size());
            #pragma omp for
            for (size_t i = 0; i < result.length(); i++) {
                size_t divider = result.length()/result.dims()[0];
                size_t index = 0;
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%result.dims()[j];
                    index += divider*indices1[j];
                    divider /= j == indices1.size() - 1 ? 1 : result.dims()[j + 1];
                }
                result(indices1) = 0;
                for (size_t j = 0; j < dimensions[dimension]; j++) {
                    std::vector<size_t> indices2 = indices1;
                    indices2[dimension] = j;
                    result(indices1) += (*this)(indices2);
                }
            }
        }
        return result;
    }

    void operator= (const Tensor<T>& src) {
        size = src.length();
        dimensions = src.dims();
        delete[] ptr;
        ptr = new T[size];
        for(size_t i = 0; i < size; i++) {
            std::vector<size_t> indices(dimensions.size());
            size_t divider = size/dimensions[0];
            size_t index = 0;
            for (size_t j = 0; j < indices.size(); j++) {
                indices[j] = (i/divider)%dimensions[j];
                index += divider*indices[j];
                divider /= j == indices.size() - 1 ? 1 : dimensions[j + 1];
            }
            ptr[index] = src(indices);
        }
    }

    void operator= (Tensor<T>&& src) {
        size = src.length();
        dimensions = src.dims();
        ptr = &src;
        src.del();
    }

    Tensor<T> operator~ () {
        std::vector<size_t> dimen = dimensions.size() == 1 ? std::vector<size_t>{size, 1} : std::vector<size_t>{dimensions[1], dimensions[0]};
        Tensor<T> result =  dimensions.size() == 1 ? Tensor<T>(vec(), dimen) : Tensor<T>(dimen);
        if (dimensions.size() > 2) {
            std::cerr << "\x1B[31mERROR: TRABSPOSE NOT SUPPORTED\033[0m\n";
            exit(1);
        }
        else if (dimensions.size() == 2) {
            #pragma omp parallel for collapse(2)
            for (size_t i = 0; i < dimensions[1]; i++) {
                for (size_t j = 0; j < dimensions[0]; j++)
                    result(i, j) = (*this)(j ,i);
            }
        }
        return result;
    }


};

//
// GREATER THAN
//
template <class T, class T1> Tensor<T> operator> (const Tensor<T>& t1, T1 a) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = t1(indices1) > a ? 1 : 0;
        }
    }
    return result;
}

template <class T, class T1> Tensor<T> operator> (T1 a, const Tensor<T>& t1) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = a > t1(indices1) ? 1 : 0;
        }
    }
    return result;
}

template <class T> Tensor<T> operator> (const Tensor<T>& t1, const Tensor<T>& t2) {
    bool first = true;
    int diff;
    if (t1.length() > t2.length()) {
        first = true;
        diff = t1.dims().size() - t2.dims().size();
        for (size_t i = 0; i < t2.dims().size(); i++) {
            if (t1.dims()[i + diff] % t2.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else if (t1.length() < t2.length()) {
        first = false;
        diff = t2.dims().size() - t1.dims().size();
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (t2.dims()[i + diff] % t1.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else {
        diff = 0;
        bool found = false;
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (!found) {
                if (t1.dims()[i] > t2.dims()[i]) {
                    found = true;
                    first = true;
                    i--;
                }
                else if (t1.dims()[i] < t2.dims()[i]) {
                    found = true;
                    first = false;
                    i--;
                }
            }
            else if (first) {
                if (t1.dims()[i] % t2.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
            else {
                if (t2.dims()[i] % t1.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
        }
    }
    Tensor<T> result =  first ? Tensor<T>(t1.dims()) : Tensor<T>(t2.dims());
    if (first) {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t2.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t1.length(); i++) {
                size_t divider = t1.length()/t1.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t1.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t2.dims()[j];
                }
                result(indices1) = t1(indices1) > t2(indices2) ? 1 : 0;
            }
        }
    }
    else {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t2.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t1.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t2.length(); i++) {
                size_t divider = t2.length()/t2.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t2.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t2.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t1.dims()[j];
                }
                result(indices1) = t1(indices2) > t2(indices1) ? 1 : 0;
            }
        }
    }
    return result;
}
//
// LESS THAN
//
template <class T, class T1> Tensor<T> operator< (const Tensor<T>& t1, T1 a) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = t1(indices1) < a ? 1 : 0;
        }
    }
    return result;
}

template <class T, class T1> Tensor<T> operator< (T1 a, const Tensor<T>& t1) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = a < t1(indices1) ? 1 : 0;
        }
    }
    return result;
}

template <class T> Tensor<T> operator< (const Tensor<T>& t1, const Tensor<T>& t2) {
    bool first = true;
    int diff;
    if (t1.length() > t2.length()) {
        first = true;
        diff = t1.dims().size() - t2.dims().size();
        for (size_t i = 0; i < t2.dims().size(); i++) {
            if (t1.dims()[i + diff] % t2.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else if (t1.length() < t2.length()) {
        first = false;
        diff = t2.dims().size() - t1.dims().size();
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (t2.dims()[i + diff] % t1.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else {
        diff = 0;
        bool found = false;
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (!found) {
                if (t1.dims()[i] > t2.dims()[i]) {
                    found = true;
                    first = true;
                    i--;
                }
                else if (t1.dims()[i] < t2.dims()[i]) {
                    found = true;
                    first = false;
                    i--;
                }
            }
            else if (first) {
                if (t1.dims()[i] % t2.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
            else {
                if (t2.dims()[i] % t1.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
        }
    }
    Tensor<T> result =  first ? Tensor<T>(t1.dims()) : Tensor<T>(t2.dims());
    if (first) {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t2.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t1.length(); i++) {
                size_t divider = t1.length()/t1.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t1.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t2.dims()[j];
                }
                result(indices1) = t1(indices1) < t2(indices2) ? 1 : 0;
            }
        }
    }
    else {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t2.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t1.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t2.length(); i++) {
                size_t divider = t2.length()/t2.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t2.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t2.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t1.dims()[j];
                }
                result(indices1) = t1(indices2) < t2(indices1) ? 1 : 0;
            }
        }
    }
    return result;
}
//
// ADDITION
//
template <class T, class T1> Tensor<T> operator+ (const Tensor<T>& t1, T1 a) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = t1(indices1) + a;
        }
    }
    return result;
}

template <class T, class T1> Tensor<T> operator+ (T1 a, const Tensor<T>& t1) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = a + t1(indices1);
        }
    }
    return result;
}

template <class T> Tensor<T> operator+ (const Tensor<T>& t1, const Tensor<T>& t2) {
    bool first = true;
    int diff;
    if (t1.length() > t2.length()) {
        first = true;
        diff = t1.dims().size() - t2.dims().size();
        for (size_t i = 0; i < t2.dims().size(); i++) {
            if (t1.dims()[i + diff] % t2.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else if (t1.length() < t2.length()) {
        first = false;
        diff = t2.dims().size() - t1.dims().size();
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (t2.dims()[i + diff] % t1.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else {
        diff = 0;
        bool found = false;
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (!found) {
                if (t1.dims()[i] > t2.dims()[i]) {
                    found = true;
                    first = true;
                    i--;
                }
                else if (t1.dims()[i] < t2.dims()[i]) {
                    found = true;
                    first = false;
                    i--;
                }
            }
            else if (first) {
                if (t1.dims()[i] % t2.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
            else {
                if (t2.dims()[i] % t1.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
        }
    }
    Tensor<T> result =  first ? Tensor<T>(t1.dims()) : Tensor<T>(t2.dims());
    if (first) {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t2.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t1.length(); i++) {
                size_t divider = t1.length()/t1.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t1.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t2.dims()[j];
                }
                result(indices1) = t1(indices1) + t2(indices2);
            }
        }
    }
    else {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t2.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t1.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t2.length(); i++) {
                size_t divider = t2.length()/t2.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t2.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t2.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t1.dims()[j];
                }
                result(indices1) = t1(indices2) + t2(indices1);
            }
        }
    }
    return result;
}
//
// SUBSTRACTION
//
template <class T, class T1> Tensor<T> operator- (const Tensor<T>& t1, T1 a) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = t1(indices1) - a;
        }
    }
    return result;
}

template <class T, class T1> Tensor<T> operator- (T1 a, const Tensor<T>& t1) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = a - t1(indices1);
        }
    }
    return result;
}

 template <class T> Tensor<T> operator- (const Tensor<T>& t1) {
        return (T)0 - t1;
    }

template <class T> Tensor<T> operator- (const Tensor<T>& t1, const Tensor<T>& t2) {
    bool first = true;
    int diff;
    if (t1.length() > t2.length()) {
        first = true;
        diff = t1.dims().size() - t2.dims().size();
        for (size_t i = 0; i < t2.dims().size(); i++) {
            if (t1.dims()[i + diff] % t2.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else if (t1.length() < t2.length()) {
        first = false;
        diff = t2.dims().size() - t1.dims().size();
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (t2.dims()[i + diff] % t1.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else {
        diff = 0;
        bool found = false;
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (!found) {
                if (t1.dims()[i] > t2.dims()[i]) {
                    found = true;
                    first = true;
                    i--;
                }
                else if (t1.dims()[i] < t2.dims()[i]) {
                    found = true;
                    first = false;
                    i--;
                }
            }
            else if (first) {
                if (t1.dims()[i] % t2.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
            else {
                if (t2.dims()[i] % t1.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
        }
    }
    Tensor<T> result =  first ? Tensor<T>(t1.dims()) : Tensor<T>(t2.dims());
    if (first) {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t2.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t1.length(); i++) {
                size_t divider = t1.length()/t1.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t1.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t2.dims()[j];
                }
                result(indices1) = t1(indices1) - t2(indices2);
            }
        }
    }
    else {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t2.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t1.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t2.length(); i++) {
                size_t divider = t2.length()/t2.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t2.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t2.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t1.dims()[j];
                }
                result(indices1) = t1(indices2) - t2(indices1);
            }
        }
    }
    return result;
}
//
// MULTIPLICATION
//
template <class T, class T1> Tensor<T> operator* (const Tensor<T>& t1, T1 a) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = t1(indices1) * a;
        }
    }
    return result;
}

template <class T, class T1> Tensor<T> operator* (T1 a, const Tensor<T>& t1) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = a * t1(indices1);
        }
    }
    return result;
}

template <class T> Tensor<T> operator* (const Tensor<T>& t1, const Tensor<T>& t2) {
    bool first = true;
    int diff;
    if (t1.length() > t2.length()) {
        first = true;
        diff = t1.dims().size() - t2.dims().size();
        for (size_t i = 0; i < t2.dims().size(); i++) {
            if (t1.dims()[i + diff] % t2.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else if (t1.length() < t2.length()) {
        first = false;
        diff = t2.dims().size() - t1.dims().size();
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (t2.dims()[i + diff] % t1.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else {
        diff = 0;
        bool found = false;
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (!found) {
                if (t1.dims()[i] > t2.dims()[i]) {
                    found = true;
                    first = true;
                    i--;
                }
                else if (t1.dims()[i] < t2.dims()[i]) {
                    found = true;
                    first = false;
                    i--;
                }
            }
            else if (first) {
                if (t1.dims()[i] % t2.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
            else {
                if (t2.dims()[i] % t1.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
        }
    }
    Tensor<T> result =  first ? Tensor<T>(t1.dims()) : Tensor<T>(t2.dims());
    if (first) {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t2.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t1.length(); i++) {
                size_t divider = t1.length()/t1.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t1.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t2.dims()[j];
                }
                result(indices1) = t1(indices1) * t2(indices2);
            }
        }
    }
    else {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t2.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t1.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t2.length(); i++) {
                size_t divider = t2.length()/t2.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t2.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t2.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t1.dims()[j];
                }
                result(indices1) = t1(indices2) * t2(indices1);
            }
        }
    }
    return result;
}
//
// DIVISION
//
template <class T, class T1> Tensor<T> operator/ (const Tensor<T>& t1, T1 a) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = t1(indices1) / a;
        }
    }
    return result;
}

template <class T, class T1> Tensor<T> operator/ (T1 a, const Tensor<T>& t1) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = a / t1(indices1);
        }
    }
    return result;
}

template <class T> Tensor<T> operator/ (const Tensor<T>& t1, const Tensor<T>& t2) {
    bool first = true;
    int diff;
    if (t1.length() > t2.length()) {
        first = true;
        diff = t1.dims().size() - t2.dims().size();
        for (size_t i = 0; i < t2.dims().size(); i++) {
            if (t1.dims()[i + diff] % t2.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else if (t1.length() < t2.length()) {
        first = false;
        diff = t2.dims().size() - t1.dims().size();
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (t2.dims()[i + diff] % t1.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else {
        diff = 0;
        bool found = false;
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (!found) {
                if (t1.dims()[i] > t2.dims()[i]) {
                    found = true;
                    first = true;
                    i--;
                }
                else if (t1.dims()[i] < t2.dims()[i]) {
                    found = true;
                    first = false;
                    i--;
                }
            }
            else if (first) {
                if (t1.dims()[i] % t2.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
            else {
                if (t2.dims()[i] % t1.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
        }
    }
    Tensor<T> result =  first ? Tensor<T>(t1.dims()) : Tensor<T>(t2.dims());
    if (first) {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t2.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t1.length(); i++) {
                size_t divider = t1.length()/t1.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t1.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t2.dims()[j];
                }
                result(indices1) = t1(indices1) / t2(indices2);
            }
        }
    }
    else {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t2.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t1.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t2.length(); i++) {
                size_t divider = t2.length()/t2.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t2.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t2.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t1.dims()[j];
                }
                result(indices1) = t1(indices2) / t2(indices1);
            }
        }
    }
    return result;
}

//
// MODULO
//
template <class T, class T1> Tensor<T> operator% (const Tensor<T>& t1, T1 a) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = t1(indices1) % a;
        }
    }
    return result;
}

template <class T, class T1> Tensor<T> operator% (T1 a, const Tensor<T>& t1) {
    Tensor<T> result = Tensor<T>(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = a % t1(indices1);
        }
    }
    return result;
}

template <class T> Tensor<T> operator% (const Tensor<T>& t1, const Tensor<T>& t2) {
    bool first = true;
    int diff;
    if (t1.length() > t2.length()) {
        first = true;
        diff = t1.dims().size() - t2.dims().size();
        for (size_t i = 0; i < t2.dims().size(); i++) {
            if (t1.dims()[i + diff] % t2.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else if (t1.length() < t2.length()) {
        first = false;
        diff = t2.dims().size() - t1.dims().size();
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (t2.dims()[i + diff] % t1.dims()[i] != 0) {
                std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                exit(1);
            }
        }
    }
    else {
        diff = 0;
        bool found = false;
        for (size_t i = 0; i < t1.dims().size(); i++) {
            if (!found) {
                if (t1.dims()[i] > t2.dims()[i]) {
                    found = true;
                    first = true;
                    i--;
                }
                else if (t1.dims()[i] < t2.dims()[i]) {
                    found = true;
                    first = false;
                    i--;
                }
            }
            else if (first) {
                if (t1.dims()[i] % t2.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
            else {
                if (t2.dims()[i] % t1.dims()[i] != 0) {
                    std::cerr << "\x1B[31mERROR: INCOMPATIBLE DIMENSIONS\033[0m\n";
                    exit(1);
                }
            }
        }
    }
    Tensor<T> result =  first ? Tensor<T>(t1.dims()) : Tensor<T>(t2.dims());
    if (first) {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t1.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t2.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t1.length(); i++) {
                size_t divider = t1.length()/t1.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t1.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t2.dims()[j];
                }
                result(indices1) = t1(indices1) % t2(indices2);
            }
        }
    }
    else {
        #pragma omp parallel
        {
            std::vector<size_t> indices1 = std::vector<size_t>(t2.dims().size());
            std::vector<size_t> indices2 = std::vector<size_t>(t1.dims().size());
            #pragma omp for
            for (size_t i = 0; i < t2.length(); i++) {
                size_t divider = t2.length()/t2.dims()[0];
                for (size_t j = 0; j < indices1.size(); j++) {
                    indices1[j] = (i/divider)%t2.dims()[j];
                    divider /= j == indices1.size() - 1 ? 1 : t2.dims()[j + 1];
                }
                for (size_t j = 0; j < indices2.size(); j++) {
                    indices2[j] = indices1[j + diff]%t1.dims()[j];
                }
                result(indices1) = t1(indices2) % t2(indices1);
            }
        }
    }
    return result;
}
//
// MATRIX MULTIPLICATION
//
template <class T> Tensor<T> mul (const Tensor<T>& t1, const Tensor<T>& t2) {
    if (t1.dims().size() != 2 || t2.dims().size() != 2) {
        std::cerr << "\x1B[31mERROR: MATRIX MULTIPLICATION NOT SUPPORTED\033[0m\n";
        exit(1);
    }
    if (t1.dims()[1] != t2.dims()[0]) {
        std::cerr << "\x1B[31mERROR: MATRIX MULTIPLICATION FAILED\033[0m\n";
        exit(1);
    }
    std::vector dimen = {t1.dims()[0], t2.dims()[1]};
    Tensor<T> result(dimen);
    #pragma omp parallel for collapse(2)
    for(size_t i = 0; i < result.dims()[0]; i++) {
        for(size_t j = 0; j < result.dims()[1]; j++) {
            result(i, j) = 0;
            for(size_t k = 0; k < t1.dims()[1]; k++) {
                result(i, j) += t1(i, k)*t2(k, j);
            }
        }
    }
    return result;
}
//
// TRIGONOMETRIC FUNCTIONS
//
template <class T> Tensor<T> cos (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = cos(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> sin (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = sin(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> tan (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = tan(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> acos (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = acos(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> asin (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = asin(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> atan (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = atan(t1(indices1));
        }
    }
    return result;
}
//
// HYPERBOLIC FUNCTIONS
//
template <class T> Tensor<T> cosh (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::cosh(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> sinh (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::sinh(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> tanh (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::tanh(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> acosh (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::acosh(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> asinh (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::asinh(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> atanh (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::atanh(t1(indices1));
        }
    }
    return result;
}
//
// EXPONENTIAL FUNCTIONS
//
template <class T> Tensor<T> exp (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::exp(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> exp2 (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::exp2(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> expm1 (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::expm1(t1(indices1));
        }
    }
    return result;
}
//
// LOGARITHMIC FUNCTIONS
//
template <class T> Tensor<T> log (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {

            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::log(t1(indices1));
            if (std::isinf(result(indices1)))result(indices1) = -std::numeric_limits<double>::max();
        }
    }
    return result;
}
template <class T> Tensor<T> log10 (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::log10(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> log2 (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::log2(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> log1p (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::log1p(t1(indices1));
        }
    }
    return result;
}
//
// ROOT FUNCTIONS
//
template <class T> Tensor<T> sqrt (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::sqrt(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> cbrt (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::cbrt(t1(indices1));
        }
    }
    return result;
}
//
// ROUNDING FUNCTIONS
//
template <class T> Tensor<T> ceil (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::ceil(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> floor (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::floor(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> round (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::round(t1(indices1));
        }
    }
    return result;
}
template <class T> Tensor<T> trunc (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::trunc(t1(indices1));
        }
    }
    return result;
}
//
// OTHER FUNCTIONS
//
template <class T, class T1> Tensor<T> min (const Tensor<T>& t1, T1 a) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::fmin(t1(indices1), a);
        }
    }
    return result;
}
template <class T, class T1> Tensor<T> max (const Tensor<T>& t1, T1 a) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::fmax(t1(indices1), a);
        }
    }
    return result;
}
template <class T> Tensor<T> abs (const Tensor<T>& t1) {
    Tensor<T> result(t1.dims());
    #pragma omp parallel
    {
        std::vector<size_t> indices1(t1.dims().size());
        #pragma omp for
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            result(indices1) = std::abs(t1(indices1));
        }
    }
    return result;
}
