#include"linearvector.h"
#include"matrix.h"
#include"QuetionA.h"
#include"evolution.h"

/*
source: The Matrix object you want to solve;
size: Number of individuals in the species (200 in the article);
max_generation: Maximum iteration times (1000 in the article);
cross_pos: Possibility of crossing (0.6 in the article);
cross_rate: Degree of crossing (0.66 in the article);
mutate_pos: Possibility of mutation (0.05 in the article);
mutate_rate: Dependence of mutation degree on the generation times (0.8 in the article);
*/
Eigenvalue::Eigenvalue(const Matrix& source, const int& size, const int& max_generation, const double& cross_pos, const double& cross_rate, const double& mutate_pos, const double& mutate_rate) {
	this->source = source;
	this->size = size;
	this->max_generation = max_generation;
	this->cross_pos = cross_pos;
	this->cross_rate = cross_rate;
	this->mutate_pos = mutate_pos;
	this->mutate_rate = mutate_rate;
}

void Eigenvalue::initiallize(std::mt19937 gen, std::uniform_real_distribution<double> u) {
	generation = 0;
	parents.clear(); children.clear();
	eigenpolyval.clear(); sum_cost.clear();
	upper_bound = 0.0;
	for (int i = 0; i < source.get_row(); i++) {
		for (int j = 0; j < source.get_column(); j++) {
			upper_bound += std::fabs(source(i,j));
		}
	}
	lower_bound = -upper_bound;
	answers.clear();
	for (int i = 0; i < size; i++) {
		parents.push_back(u(gen)*(upper_bound-lower_bound)+lower_bound);
	}
}
double Eigenvalue::solve_eigenpolyval(double lam) {
	return std::fabs((source-lam*identity_matrix(source.get_column())).determinant());
}
void Eigenvalue::solve_total() {
	double max_eigenpolyval = 0.0;
	eigenpolyval.clear(); sum_cost.clear();
	for (int i = 0; i < size; i++) {
		eigenpolyval.push_back(solve_eigenpolyval(parents[i]));
		if (eigenpolyval[i] > max_eigenpolyval) max_eigenpolyval = eigenpolyval[i];
	}
	sum_cost.push_back(max_eigenpolyval - eigenpolyval[0]);
	for (int i = 1; i < size; i++) {
		sum_cost.push_back(sum_cost[i-1] + max_eigenpolyval - eigenpolyval[i]);
	}
}
int Eigenvalue::choose(double proportion) {
	for (int i = 0; i < size; i++) {
		if (sum_cost[i] > proportion * sum_cost[size - 1]) return i;
	}
}
void Eigenvalue::mutate(int idx, bool opt, double ran) {
	double tmp = children[idx];
	if (opt) {
		children[idx] += (upper_bound - tmp) * (1.0 - std::pow(ran, mutate_rate * (double)generation / max_generation));
	}
	else {
		children[idx] -= (tmp - lower_bound) * (1.0 - std::pow(ran, mutate_rate * (double)generation / max_generation));
	}
}
void Eigenvalue::crossover(int ida, int idb, bool opt) {
	if (opt) {
		children.push_back(cross_rate * parents[idb] + (1.0 - cross_rate) * parents[ida]);
		children.push_back(cross_rate * parents[ida] + (1.0 - cross_rate) * parents[idb]);
	}
	else {
		children.push_back(parents[ida]);
		children.push_back(parents[idb]);
	}
}
void Eigenvalue::produce_next(std::mt19937 gen, std::uniform_real_distribution<double> u) {
	children.clear();
	generation++;
	solve_total();
	int ida, idb;
	for (int i = 0; i < (size >> 1); i++) {
		ida = choose(u(gen));
		idb = choose(u(gen));
		crossover(ida, idb, u(gen)<cross_pos);
	}
	for (int i = 0; i < size; i++) if (u(gen) < mutate_pos) {
		mutate(i,u(gen)<0.5,u(gen));
	}
	std::swap(parents, children);
}
void Eigenvalue::solve_eigenvalue() {
	std::mt19937 gen(time(NULL));
	std::uniform_real_distribution<double> u(0.0, 1.0);
	if (source.get_column() != source.get_row()) {
		std::cerr << "[ERROR]Non-square matrices have no eigenvalue!!" << std::endl;
		return;
	}
	initiallize(gen,u);
	for (int i = 0; i < max_generation; i++) {
		produce_next(gen, u);
	}
	std::sort(parents.begin(), parents.end());

	for (int i = 0; i < size; i++) {
		if (solve_eigenpolyval(parents[i]) < eps2&&(answers.empty()||parents[i]-*(answers.end()-1)>eps))	answers.push_back(parents[i]);
	}
}
