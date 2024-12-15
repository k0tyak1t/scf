#ifndef MO_H
#define MO_H
#include "matrix.h"

class MO {
public:
    MO(const int);
    MO();
    ~MO();
    matrix C; // MO to AO expansion. MO in columns
    const double get_mo_energy(const int) const;
    const int get_size() const;
    int set_mo_energy(const int, const double);
    int set_total_energy(const double);
    const double get_total_energy() const;
    int get_irrep_characters(const int&) const;
    bool set_c2v(int* symmAO, const double& limit);
    void set_mo_energies(const double*);

private:
    int n; // basis dimentionality
    double* mo_energies;
    double total_energy;
    int* irrep;
};

#endif // MO_H
