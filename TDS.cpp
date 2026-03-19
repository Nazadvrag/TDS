#include "TDS.h"

TDS::TDS()
{
	size = 1;
	method = nullptr;
	a = new double[size];
	b = new double[size];
	c = new double[size];
	f = new double[size];
}

TDS::TDS(double* _a, double* _b, double* _c, double* _f, int n)
{
	size = n;
	method = nullptr;
	a = new double[size];
	b = new double[size];
	c = new double[size];
	f = new double[size];
	for (int i = 0; i < n; i++)
	{
		a[i] = _a[i];
		b[i] = _b[i];
		c[i] = _c[i];
		f[i] = _f[i];
	}
}

TDS::~TDS()
{
	delete[] f;
	delete[] c;
	delete[] b;
	delete[] a;
	delete method;
}

TDS::TDS(const TDS& other) :
	size(other.size),
	method(nullptr),
	a(nullptr),
	b(nullptr),
	c(nullptr),
	f(nullptr)
{
	if (size > 0) {
		a = new double[size];
		b = new double[size];
		c = new double[size];
		f = new double[size];

		for (int i = 0; i < size; i++) {
			a[i] = other.a[i];
			b[i] = other.b[i];
			c[i] = other.c[i];
			f[i] = other.f[i];
		}
	}
}

TDS& TDS::operator=(const TDS& other) {
	if (this == &other) return *this;

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] f;
	delete method;

	size = other.size;
	method = nullptr;

	if (size > 0) {
		a = new double[size];
		b = new double[size];
		c = new double[size];
		f = new double[size];
		
		for (int i = 0; i < size; ++i) 
		{
			a[i] = other.a[i];
			b[i] = other.b[i];
			c[i] = other.c[i];
			f[i] = other.f[i];
		}
	 }
	else {
		a = b = c = f = nullptr;
	}
	
	return *this;
}

void TDS::setSize(int newSize)
{
	size = newSize;

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] f;

	a = new double[size];
	b = new double[size];
	c = new double[size];
	f = new double[size];
}

double* TDS::randomFill(int range)
{
	double* x = new double[size];
	a[0] = 0;
	c[size - 1] = 0;
	for (int i = 0; i < size; i++)
	{
		x[i] = static_cast<double>(rand()) / RAND_MAX * (range * 2) - range;
		if (i!=0) a[i] = static_cast<double>(rand()) / RAND_MAX * (range * 2) - range;
		b[i] = static_cast<double>(rand()) / RAND_MAX * (range * 2) - range;
		if (i != size-1) c[i] = static_cast<double>(rand()) / RAND_MAX * (range * 2) - range;
	}
	f[0] = b[0] * x[0] + c[0] * x[1];
	f[size - 1] = a[size - 1] * x[size - 2] + b[size - 1] * x[size - 1];
	for (int i = 1; i < size-1; i++)
	{
		f[i] = a[i] * x[i - 1] + b[i] * x[i] + c[i] * x[i + 1];
	}
	return x;
	delete[] x;
}

double* TDS::testFill()
{
	a[0] = 0;
	c[size - 1] = 0;
	double* x = new double[size];
	for (int i = 0; i < size; i++)
	{
		if (i != 0) a[i] = 1;
		if (i != size - 1) c[i] = 1;
		b[i] = 100;
		f[i] = a[i] + c[i] + b[i];
		x[i] = 1;
	}
	return x;
}
double* TDS::testFill2()
{
	a[0] = 0;
	c[size - 1] = 0;
	double* x = new double[size];
	for (int i = 0; i < size; i++)
	{
		if (i != 0) a[i] = 1;
		if (i != size - 1) c[i] = 1;
		b[i] = i%2==1?2:3;
		b[0] = 2;
		b[size - 1] = 3;
		x[i] = i % 2 == 0 ? 2 : 1;
	}
	f[0] = b[0] * x[0] + c[0] * x[1];
	f[size - 1] = a[size - 1] * x[size - 2] + b[size - 1] * x[size - 1];
	for (int i = 1; i < size - 1; i++) f[i] = a[i] * x[i - 1] + b[i] * x[i] + c[i] * x[i + 1];
	return x;
}
int TDS::getSize()
{
	return size;
}

void TDS::setMethod(SolvingMethod* newMethod)
{
	method = move(newMethod);
}

double ThomasMethod::solve(double* a, double* b, double* c, double* f, int n, int num_threads, double*& x)
{
	double start = omp_get_wtime();
	for (int i = 0; i < n-1; i++)
	{
		c[i] /= b[i];
		f[i] /= b[i];
		b[i] = 1;
		double coef = a[i + 1];
		a[i + 1] = 0;
		b[i + 1] = b[i+1] - coef * c[i];
		f[i + 1] = f[i+1] - coef * f[i];
	}
	//double* x = new double[n];
	x[n - 1] = f[n - 1] / b[n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = (f[i] - c[i] * x[i + 1]) / b[i];
	}
	double time = omp_get_wtime() - start;
	return time;
}

double TDS::solve(int num_threads, double*& x)
{
	double result;
	if (method) 
	{
		result = method->solve(a, b, c, f, size, num_threads, x);	
	}
	else 
	{
		std::cout << "Метод не установлен!" << endl;
	}
	return result;
}

double PartitionMethod::solve(double* a, double* b, double* c, double* f, int size, int num_threads, double*& x)
{
	omp_set_num_threads(num_threads);
	double t_start = omp_get_wtime();
	int* parts = new int[num_threads + 1];
	parts[0] = 0;
	parts[num_threads] = size;
	for (int i = 1; i < num_threads; i++)
	{
		parts[i] = size / num_threads * i;
	}
#pragma omp parallel
	{
		int thread_num = omp_get_thread_num();
		if (thread_num == 0)
		{
			for (int i = parts[0]; i < parts[1]-1; i++)
			{
				c[i] /= b[i];
				f[i] /= b[i];
				b[i] = 1;
				double coef = a[i + 1];
				a[i + 1] = 0;
				b[i + 1] = b[i + 1] - coef * c[i];
				f[i + 1] = f[i + 1] - coef * f[i];
			}
		}
		else
		{
			for (int i = parts[thread_num]; i < parts[thread_num + 1]-1; i++)
			{
				a[i] /= b[i];
				c[i] /= b[i];
				f[i] /= b[i];
				b[i] = 1;
				double coef = a[i + 1];
				a[i + 1] = -a[i] * coef;
				b[i + 1] = b[i + 1] - coef * c[i];
				f[i + 1] = f[i + 1] - coef * f[i];
			}	
		}
		for (int i = parts[thread_num + 1] - 2; i > parts[thread_num]; i--)
		{
			double coef = c[i - 1];
			c[i - 1] = -coef * c[i];
			a[i - 1] -= coef * a[i];
			f[i - 1] -= coef * f[i];
		}
	}
	for (int i = 1; i < num_threads; i++)
	{
		double coef = c[parts[i]-1];
		c[parts[i] - 1] = -coef * c[parts[i]];
		b[parts[i] - 1] -= coef * a[parts[i]];
		f[parts[i] - 1] -= coef * f[parts[i]];
	}
	double* sa = new double[num_threads];
	double* sb = new double[num_threads];
	double* sc = new double[num_threads];
	double* sf = new double[num_threads];
	sa[0] = 0;
	sc[num_threads - 1] = 0;
	for (int i = 0; i < num_threads; i++)
	{
		if (i != 0) sa[i] = a[parts[i + 1] - 1];
		sb[i] = b[parts[i + 1] - 1];
		if (i!=num_threads-1) sc[i] = c[parts[i + 1] - 1];
		sf[i] = f[parts[i + 1] - 1];
	}
	TDS* sm = new TDS(sa, sb, sc, sf, num_threads);
	sm->setMethod(new ThomasMethod());
	double* sx = new double[num_threads];
	sm->solve(1, sx);
	delete sm;
	delete[] sa;
	delete[] sb;
	delete[] sc;
	delete[] sf;
	for (int i = 0; i < num_threads; i++)
	{
		x[parts[i + 1] - 1] = sx[i];
	}
	delete[] sx;
	int count=0;
#pragma omp parallel
	{
		int thread_num = omp_get_thread_num();
		if (thread_num == 0)
		{
			for (int i = parts[1] - 2; i >= parts[0];i--)
			{
				x[i] = (f[i] - c[i] * x[parts[1]-1]) / b[i];
			}
		}
		else
		{
			for (int i = parts[thread_num+1] - 2; i > parts[thread_num]-1; i--)
			{
				x[i] = (f[i] - c[i] * x[parts[thread_num+1]-1] - a[i]* x[parts[thread_num] - 1]) / b[i];
			}
		}
	}
	count = 0;
	double time = omp_get_wtime() - t_start;
	return time;
	delete[] parts;
}
/*
double* reductionMethod::solve(double* a, double* b, double* c, double* f, int n)
{
	omp_set_num_threads(12);
	double* x;
	if (n < 5) x = ThomasMethod().solve(a, b, c, f, n);
	else
	{
		int new_n = n/2 + n%2;
		double* new_a = new double[new_n];
		double* new_b = new double[new_n];
		double* new_c = new double[new_n];
		std::cout << "a" << endl;
		double* new_f = new double[new_n];
		a[0] = 0;
		c[new_n - 1] = 0;
# pragma omp parallel for
		for (int i = 0; i < new_n; i++)
		{
			if (i != 0 && i!=new_n-1)
			{
				double up_coef = -a[2 * i] / b[2 * i - 1];
				double down_coef = -c[2 * i] / b[2 * i + 1];
				new_a[i] = up_coef * a[2 * i - 1];
				new_c[i] = down_coef * c[2 * i + 1];
				new_b[i] = b[2 * i] + up_coef * c[2 * i - 1] + down_coef * a[2 * i + 1];
				new_f[i] = f[2 * i] + up_coef * f[2 * i - 1] + down_coef * f[2 * i + 1];
			}
			else if (i == 0)
			{
				double down_coef= -c[2 * i] / b[2 * i + 1];
				new_c[i] = down_coef * c[2 * i + 1];
				new_b[i] = b[2 * i] + down_coef * a[2 * i + 1];
				new_f[i] = f[2 * i] + down_coef * f[2 * i + 1];
			}
			else if (new_n % 2 == 1)
			{
				double up_coef = -a[2 * i] / b[2 * i - 1];
				new_a[i] = up_coef * a[2 * i - 1];
				new_b[i] = b[2 * i] + up_coef * c[2 * i - 1];
				new_f[i] = f[2 * i] + up_coef * f[2 * i - 1];
			}
			else
			{
				double up_coef = -a[2 * i] / b[2 * i - 1];
				double down_coef = -c[2 * i] / b[2 * i + 1];
				new_a[i] = up_coef * a[2 * i - 1];
				new_b[i] = b[2 * i] + up_coef * c[2 * i - 1] + down_coef * a[2 * i + 1];
				new_f[i] = f[2 * i] + up_coef * f[2 * i - 1] + down_coef * f[2 * i + 1];
			}
		}
		double* new_x = reductionMethod().solve(new_a, new_b, new_c, new_f, new_n);
		x = new double[n];
		for (int i = 0; i < n; i++)
		{
			x[2 * i] = new_x[i];
		}
# pragma omp parallel for
		for (int i = 1; i < n; i += 2)
		{
			if (n % 2 == 0 && i == n - 1) x[i] = (f[i] - a[i] * x[i - 1]) / b[i];
			else x[i] = (f[i] - c[i] * x[i + 1] - a[i] * x[i - 1])/b[i];
		}
	}
	return x;
}*/
double ReductionMethod::solve(double* a, double* b, double* c, double* f, int n, int num_threads, double*& x)
{
	omp_set_num_threads(num_threads);
	double t_start = omp_get_wtime();
	int nlog = log2(n + 1);
	for (int q = 1; q < nlog-1; q++)
	{
		int step = pow(2, q);
		int gap = step / 2;
#pragma omp parallel for
			for (int i = step - 1; i < n; i += step)
			{
				if (i == step - 1)
				{
					double coef = a[i] / b[i - gap];
					a[i] = 0;
					b[i] -= coef * c[i - gap];
					f[i] -= coef * f[i - gap];
					coef = c[i] / b[i + gap];
					c[i] = -coef * c[i + gap];
					b[i] -= coef * a[i + gap];
					f[i] -= coef * f[i + gap];
				}
				else
					if (i == n - step + 1)
					{
						double coef = a[i] / b[i - gap];
						a[i] = -coef * a[i - gap];
						b[i] -= coef * c[i - gap];
						f[i] -= coef * f[i - gap];
						coef = c[i] / b[i + gap];
						c[i] = 0;
						b[i] -= coef * a[i + gap];
						f[i] -= coef * f[i + gap];
					}
					else
					{
						double coef = a[i] / b[i - gap];
						a[i] = -coef * a[i - gap];
						b[i] -= coef * c[i - gap];
						f[i] -= coef * f[i - gap];
						coef = c[i] / b[i + gap];
						c[i] = -coef * c[i + gap];
						b[i] -= coef * a[i + gap];
						f[i] -= coef * f[i + gap];
					}
			}
	}
	int mid = n/2;
	int gap = mid / 2 + 1;
	double coef = a[mid] / b[mid - gap];
	a[mid] = 0;
	b[mid] -= coef * c[mid - gap];
	f[mid] -= coef * f[mid - gap];
	coef = c[mid] / b[mid + gap];
	c[mid] = 0;
	b[mid] -= coef * a[mid + gap];
	f[mid] -= coef * f[mid + gap];
	x[mid] = f[mid] / b[mid];
	for (int q = nlog - 1; q > 0; q--)
	{
		int step = pow(2, q);
		int gap = step / 2;
#pragma omp parallel for
			for (int i = gap - 1; i < n; i += step)
			{
				if (i == gap - 1) x[i] = (f[i] - c[i] * x[i + gap]) / b[i];
				else if (i == n - gap + 1) x[i] = (f[i] - a[i] * x[i - gap]) / b[i];
				else x[i] = (f[i] - a[i] * x[i - gap] - c[i] * x[i + gap]) / b[i];
			}
	}
	double time = omp_get_wtime() - t_start;
	return time;
}