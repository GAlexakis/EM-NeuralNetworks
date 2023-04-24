#include <iostream>
#include <vector>

template<class T> class Tensor {
private:
    size_t size;
    std::vector<size_t> dimensions;
    T* ptr;

public:
    Tensor () {
        dimensions = {};
        ptr = new T{};
    }

    Tensor (std::vector<size_t> dims) {
        size = 1;
        dimensions = dims;
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
        std::cout << "copy cosntructor\n";
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
        std::cout << "move cosntructor\n";
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

    size_t operator[] (size_t dimension) const {
        if (dimension >= dimensions.size() || dimension < 0) {
            std::cerr << "\x1B[31mERROR: OUT OF DIMENSION BOUNDS\033[0m\n";
            exit(1);
        }
        return dimensions[dimension];
    }

    void operator= (const Tensor<T>& src) {
        std::cout << "copy assignment\n";
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
        std::cout << "move assignment\n";
        size = src.length();
        dimensions = src.dims();
        ptr = &src;
        src.del();
    }
};

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
    std::vector<size_t> indices1 = first ? std::vector<size_t>(t1.dims().size()) : std::vector<size_t>(t2.dims().size());
    std::vector<size_t> indices2 = first ? std::vector<size_t>(t2.dims().size()) : std::vector<size_t>(t1.dims().size());
    if (first) {
        for (size_t i = 0; i < t1.length(); i++) {
            size_t divider = t1.length()/t1.dims()[0];
            size_t index = 0;
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t1.dims()[j];
                index += divider*indices1[j];
                divider /= j == indices1.size() - 1 ? 1 : t1.dims()[j + 1];
            }
            for (size_t j = 0; j < indices2.size(); j++) {
                indices2[j] = indices1[j + diff]%t2.dims()[j];
            }
            result(indices1) = t1(indices1) + t2(indices2);
        }
    }
    else {
        for (size_t i = 0; i < t2.length(); i++) {
            size_t divider = t2.length()/t2.dims()[0];
            size_t index = 0;
            for (size_t j = 0; j < indices1.size(); j++) {
                indices1[j] = (i/divider)%t2.dims()[j];
                index += divider*indices1[j];
                divider /= j == indices1.size() - 1 ? 1 : t2.dims()[j + 1];
            }
            for (size_t j = 0; j < indices2.size(); j++) {
                indices2[j] = indices1[j + diff]%t1.dims()[j];
            }
            result(indices1) = t1(indices2) + t2(indices1);
        }
    }
    return result;
}

int main(int argc, char** argv) {
    std::vector<int> values1 = {1,2,3,4};
    std::vector<int> values2 = {1,2,3,4,5,6,7,8};
    std::vector<size_t> dims1 = {2,2,1};
    std::vector<size_t> dims2 = {2,2,2};
    Tensor<int> t1(values1, dims1);
    Tensor<int> t2(values2, dims2);
    Tensor<int> t3 = t1 + t2;

    t1.print();
    t2.print();
    t3.print();
}