#include<iostream>
#include<fstream>
#include<Windows.h>
#include<omp.h>
#include<cmath>
#include<time.h>
using namespace std;
class SolvingMethod
{
public:
	virtual ~SolvingMethod() {}
	virtual double solve(double* a, double* b, double* c, double* f, int n, int num_threads, double*& x) = 0;
};
class ThomasMethod:public SolvingMethod
{
public:
	double solve(double* a, double* b, double* c, double* f, int n, int num_threads, double*& x) override;
};
class PartitionMethod : public SolvingMethod
{
public:
	double solve(double* a, double* b, double* c, double* f, int n, int num_threads, double*& x) override;
};
class reductionMethod : public SolvingMethod
{
public:
	double solve(double* a, double* b, double* c, double* f, int n, int num_threads, double*& x) override;
};
class ReductionMethod : public SolvingMethod
{
public:
	double solve(double* a, double* b, double* c, double* f, int n, int num_threads, double*& x) override;
};
class TDS
{
	double* a;
	double* b;
	double* c;
	double* f;
	int size;
	SolvingMethod* method;
public:
	TDS();
	TDS(double* _a, double* _b, double* _c, double* f, int n);
	TDS(const TDS& other);
	~TDS();
	TDS& operator=(const TDS& other);
	void setSize(int newSize);
	double* randomFill(int range);
	double* testFill();
	double* testFill2();
	int getSize();
	void setMethod(SolvingMethod* newMethod);
	double solve(int num_threads, double*& x);
};

