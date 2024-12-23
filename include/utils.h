#ifndef UTILS_H
#define UTILS_H
#include "mo.h"
#include <deque>
#include <iostream>
#include <string>
#include <vector>

void print_orbitals(const MO&);
void print_energies(const MO&);
void display_progress(int, const std::string&);
double dot_product(const std::vector<double>&, const std::vector<double>&);

#endif // UTILS_H
