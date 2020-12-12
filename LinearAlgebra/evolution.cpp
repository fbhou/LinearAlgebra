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

void initiallize(Eigenvalue e) {
	e.generation = 0;
	e.parents.clear(); e.children.clear();
	e.eigenpolyval.clear(); e.sum_cost.clear();
	e.upper_bound = 0.0;
	for (int i = 0; i < e.source.get_row(); i++) {
		for (int j = 0; j < e.source.get_column(); j++) {
			e.upper_bound += std::fabs(e.source(i,j));
		}
	}
	e.lower_bound = -e.upper_bound;
	e.answers.clear();
}
double Eigenvalue::solve_eigenpolyval(double lam) {
	tmp = source;
	for (int i = 0; i < source.get_row(); i++) {
		tmp(i, i) -= lam;
	}
	return std::fabs(tmp.determinant());
}
void solve_total(Eigenvalue e) {
	double max_eigenpolyval = 0.0;
	e.eigenpolyval.clear(); e.sum_cost.clear();
	for (int i = 0; i < e.size; i++) {
		e.eigenpolyval.push_back(e.solve_eigenpolyval(e.parents[i]));
		if (e.eigenpolyval[i] > max_eigenpolyval) max_eigenpolyval = e.eigenpolyval[i];
	}
	e.sum_cost.push_back(e.eigenpolyval[0]);
	for (int i = 1; i < e.size; i++) {
		e.sum_cost.push_back(max_eigenpolyval-e.eigenpolyval[i]);
	}
}
int Eigenvalue::choose(double proportion) {
	for (int i = 0; i < size; i++) if (sum_cost[i] > proportion * sum_cost[size - 1]) return i;
}
void mutate(Eigenvalue e, int idx, bool opt, double ran) {
	double tmp = e.children[idx];
	if (opt) {
		e.children[idx] += (e.upper_bound - tmp) * (1.0 - std::pow(ran, e.mutate_rate * (double)e.generation / e.max_generation));
	}
	else {
		e.children[idx] += (tmp - e.lower_bound) * (1.0 - std::pow(ran, e.mutate_rate * (double)e.generation / e.max_generation));
	}
}
void crossover(Eigenvalue e, int ida, int idb, bool opt) {
	if (opt) {
		e.children.push_back(e.cross_rate * e.parents[idb] + (1.0 - e.cross_rate) * e.parents[ida]);
		e.children.push_back(e.cross_rate * e.parents[ida] + (1.0 - e.cross_rate) * e.parents[idb]);
	}
	else {
		e.children.push_back(e.parents[ida]);
		e.children.push_back(e.parents[idb]);
	}
}
void produce_next(Eigenvalue e, std::mt19937 gen, std::uniform_real_distribution<double> u) {
	e.children.clear();
	e.generation++;
	solve_total(e);
	int ida, idb;
	for (int i = 0; i < (e.size >> 1); i++) {
		ida = e.choose(u(gen));
		idb = e.choose(u(gen));
		crossover(e, ida, idb, u(gen)<e.cross_pos);
	}
	for (int i = 0; i < e.size; i++) if (u(gen) < e.mutate_pos) {
		mutate(e,i,u(gen)<0.5,u(gen));
	}
	std::swap(e.parents, e.children);
}
void Eigenvalue::solve_eigenvalue() {
	std::mt19937 gen(time(NULL));
	std::uniform_real_distribution<double> u(0.0, 1.0);
	if (source.get_column() != source.get_row()) {
		std::cerr << "Baka!! Non-square matrices have no eigenvalue!!" << std::endl;
		return;
	}
	initiallize(*this);
	std::sort(parents.begin(), parents.end());
	answers.push_back(parents[0]);
	for (int i = 1; i < size; i++) if(parents[i]-*answers.end() > eps){
		answers.push_back(parents[i]);
	}
}