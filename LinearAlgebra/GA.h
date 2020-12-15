#pragma once
#ifndef EVOLUTION_H
#define EVOLUTION_H
#include"linearvector.h"
#include"matrix.h"
#include"QuetionA.h"
#include<iostream>
#include<cmath>
#include<algorithm>
#include<vector>
#include<random>
#include<ctime>

class Eigenvalue {
private:
	Matrix source;
	Matrix tmp;
	int size,generation,max_generation;
	std::vector<double> parents;
	std::vector<double> children;
	std::vector<double> eigenpolyval;
	std::vector<double> sum_cost;
	double cross_pos, cross_rate, mutate_pos, mutate_rate;
	double eps = 1e-3, eps2=1e-2;
public:

	double upper_bound, lower_bound;
	std::vector<double> answers;
	
	Eigenvalue(const Matrix&, const int&, const int&, const double&, const double&, const double&, const double&);
	~Eigenvalue()=default;

	void initiallize(std::mt19937, std::uniform_real_distribution<double>);
	double solve_eigenpolyval(double x);
	void solve_total();
	int choose(double);
	void mutate(int, bool, double);
	void crossover(int, int, bool);
	void produce_next(std::mt19937, std::uniform_real_distribution<double>);
	void solve_eigenvalue();
};


#endif

