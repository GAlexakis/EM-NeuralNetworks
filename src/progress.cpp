#include "progress.hpp"

void restore_terminal (int signum) {
    std::cout << "\x1b[0m" << "\x1b[?25h" << std::endl;
    exit(signum);
}

Progress::Progress (std::vector<int> length) {
    this->length = length;
    pos.reserve(length.size());
    bar.reserve(length.size());
    for (int i = 0; i < length.size(); i++) {
        pos.emplace_back(0);
        bar.emplace_back("");
    }
    std::cout << "\x1b[?25l";
}

Progress::~Progress () {
    std::cout << "\x1b[" << length.size() << 'E' << "\x1b[0m" << "\x1b[?25h" << std::flush;
}
void Progress::show (std::vector<std::string> msg) {
    std::cout << std::flush;
    for (int i = 0; i < length.size(); i++) {
        if (pos[i] < length[i])
            std::cout << "\x1b[33m"  << '[' << std::setw(length[i]) << std::left << bar[i] + '>' << ']' << msg[i] << '\n';
        else
            std::cout << "\x1b[32m"  << '[' << std::setw(length[i]) << std::left << std::string(length[i], '=') << ']' << msg[i] << '\n';
    }
    std::cout << "\x1b[" << length.size() << 'F';
}

void Progress::update (int index) {
    pos[index]++;
    bar[index]+= '=';
}

void Progress::reset (int index) {
    pos[index] = 0;
    bar[index] = "";
    std::cout << std::flush;
}