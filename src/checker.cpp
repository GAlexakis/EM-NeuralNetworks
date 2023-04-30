#include "checker.hpp"

bool is_number(const std::string& s) {
    try {
        std::stod(s);
    }
    catch(...) {
        return false;
    }
    return true;
}