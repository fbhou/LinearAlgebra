#include<iostream>
#include"QuetionA.h"
#include"QuetionB.h"
#include"evolution.h"
int main()
{
//	std::cout<<"////////////////////////////////////////QUETION1////////////////////////////////////////"<<std::endl;
//	QuestionA::main();
//	std::cout<<std::endl<<std::endl;
//	std::cout<<"////////////////////////////////////////QUETION2////////////////////////////////////////"<<std::endl;
//	QuestionB::main();
	double a[] = { 1,0,0,0,2,0,0,0,3 };
	Matrix A = Matrix((unsigned int)3, (unsigned int)3, a);
	Eigenvalue eigenv = Eigenvalue(A, 500, 10000, 0.6, 0.66, 0.05, 0.8);
	eigenv.solve_eigenvalue();
	printf("%lld\n", eigenv.answers.size());
	for (int i = 0; i < eigenv.answers.size(); i++) {
		printf("%.4lf\n",eigenv.answers[i]);
	}
	return 0;
}
