#ifndef RHF_H
#define RHF_H
#include <vector>
#include <deque>
#include "matrix.h"
#include "standard_matrices.h"
#include "mo.h"

class RHF {
public:
  RHF(standard_matrices&, MOs&);
	~RHF();
	bool get_convergency();
	void validate_convergency(const int);
  void update_energy();
  matrix transform_matrix(const matrix&); // todo
  matrix transform_exp_coefs(); // todo 
  void calc_density();
  void transform_density();
  void calc_eri_matrix();
  void calc_fock();
  void calc_exp_coefs();
  void print_iteration(const int);
  void core_guess();
  void direct_iteration();
	void roothan_hall();

private:
	double etol;
	int max_iter, diis_size;
	bool is_converged;
	standard_matrices& std_m;
	MOs& mo;
  double E_old, E_new;
  double* evec;
  double* eval;
  matrix density, eri_matrix, fock_matrix;
  std::deque<matrix> error_buffer;
  std::vector<double> diis_coefs;
};

#endif
