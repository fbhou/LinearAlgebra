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
	double eps = 1e-6;
public:

	double upper_bound, lower_bound;
	std::vector<double> answers;
	
	Eigenvalue(const Matrix&, const int&, const int&, const double&, const double&, const double&, const double&);
	~Eigenvalue()=default;

	friend void initiallize(Eigenvalue);
	double solve_eigenpolyval(double x);
	friend void solve_total(Eigenvalue);
	int choose(double);
	friend void mutate(Eigenvalue, int, bool, double);
	friend void crossover(Eigenvalue, int, int, bool);
	friend void produce_next(Eigenvalue, std::mt19937, std::uniform_real_distribution<double>);
	void solve_eigenvalue();
};


#endif

