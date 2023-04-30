#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <signal.h>

void restore_terminal (int signum);

class Progress {
private:
    std::vector<int> length;
    std::vector<int> pos;
    std::vector<std::string> bar;
public:
    Progress (std::vector<int> length);
    ~Progress ();
    void show (std::vector<std::string> msg);
    void update (int index);
    void reset (int index);
};