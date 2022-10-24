#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <stdio.h>

using namespace std;

double e = 1E-6;
double e1 = 1E-3;

using namespace std;

void Jacobi(double** Ak, double** X, int m);
void Norm2(double* X, int n);
void MxM(double** A, double** B, double** X, int n);
void MxN(double** A, double* Y, double* Y_new, int n);
double Apos_est(double** A, double* Y, double lambda, int n);
void Skalar(double** A, double* sol, int n);
void Degree(double** A, double* sol, int n);
void rev_Skalar(double** A, double* X, double lambda, int n);
void Gauss(double** A1, double* X, int m, int n);
void Vil(double** A, double* sol, int n);
double Eytken(double a, double b, double c);


int main()
{
	setlocale(LC_ALL, "Russian");

	cout << "Задание №16" << endl;
	cout << "1) Проблема собственных значений" << endl;
	cout << "2) Частичная проблема собственных значений" << endl;
	cout << "Вариант 2" << endl << endl;

	int m = 3;// n = 4;
	double* solution = new double[m];

	double** A = new double* [m]; // создаём массивы
	for (int i = 0; i < m; i++)
		A[i] = new double[m];

	A[0][0] = -0.82005;
	A[0][1] = -0.13542;
	A[0][2] = 0.26948;
	A[1][0] = -0.13542;
	A[1][1] = 0.51486;
	A[1][2] = 0.02706;
	A[2][0] = 0.26948;
	A[2][1] = 0.02706;
	A[2][2] = -0.83365;

	//A[0][0] = 6.0;
	//A[0][1] = 4.0;
	//A[0][2] = -2.0;
	//A[1][0] = -2.0;
	//A[1][1] = 0.0;
	//A[1][2] = 2.0;
	//A[2][0] = 6.0;
	//A[2][1] = 6.0;
	//A[2][2] = -2.0;

	cout << "Матрица А:" << endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	

	double** Ak = new double* [m]; // создаём массивы
	for (int i = 0; i < m; i++)
		Ak[i] = new double[m];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			Ak[i][j] = A[i][j];
		}
	}

	double** X = new double* [m]; // создаём массивы
	for (int i = 0; i < m; i++)
		X[i] = new double[m] {0};
	for (int i = 0; i < m; i++)
		X[i][i] = 1;

	// Находим собственные числа методом Якоби
	Jacobi(Ak, X, m);

	cout << "1) Найдём методом Якоби все собственные числа и собственные вектора матрицы А: " << endl;
	cout << "Матрица A(k):" << endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cout << setprecision(16) << Ak[i][j] << " ";
		}
		cout << endl;
	}

	for (int i = 0; i < m; i++)
		solution[i] = X[i][0];
	Norm2(solution, m);
	cout << endl<< "Вектор 1 сч" << endl;
	for (int i = 0; i < m; i++)
	{
		cout << setprecision(16) << solution[i] << " ";
	}

	for (int i = 0; i < m; i++)
		solution[i] = X[i][1];
	Norm2(solution, m);
	cout << endl<< "Вектор 2 сч" << endl;
	for (int i = 0; i < m; i++)
	{
		cout << setprecision(16) << solution[i] << " ";
	}

	for (int i = 0; i < m; i++)
		solution[i] = X[i][2];
	Norm2(solution, m);
	cout << endl<<"Вектор 3 сч" << endl;
	for (int i = 0; i < m; i++)
	{
		cout << setprecision(16) << solution[i] << " ";
	}

	cout << endl << "Матрица X:" << endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cout << setprecision(16) << X[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	
	double* sol_deg = new double[m + 2];
	double* sol_skal = new double[m + 2];
	double* sol_rev_skal = new double[m + 2];
	double* sol_vil = new double[100];

	// найдем степенным методом максимальное по модулю собственное число и его вектор
	
	Degree(A, sol_deg, m);
	for (int i = 0; i < m; i++)
	{
		solution[i]=sol_deg[i+2];
	}
	Norm2(solution, m);

	cout << endl<< endl;
	cout << "2) Максимальное по модулю собственное число, найденное степенным методом: " << sol_deg[1] << endl;
	cout << "Соответствующий ему собственный вектор: " << endl;
	for (int i = 0; i < m; i++)
	{
		cout << setprecision(16) << solution[i] << " ";
	}
	cout << endl << endl;

	// найдем методом скалярных произведений максимальное по модулю собственное число и его вектор

	Skalar(A, sol_skal, m);

	for (int i = 0; i < m; i++)
	{
		solution[i] = sol_skal[i + 2];
	}
	Norm2(solution, m);

	cout << endl;
	cout << "3) Максимальное по модулю собственное число, найденное методом скалярных произведений: " << sol_skal[1] << endl;
	cout << "Соответствующий ему собственный вектор: " << endl;
	for (int i = 0; i < m; i++)
	{
		cout << setprecision(16)<< solution[i] << " ";
	}
	cout << endl << endl;

	// найдем противоположную границу спектра собственных чисел и вектор

	rev_Skalar(A, sol_rev_skal, sol_skal[1], m);
	//for (int i = 0; i < m + 2; i++)
	//{
	//	cout << setprecision(16) << sol_rev_skal[i] << " ";
	//}
	for (int i = 0; i < m; i++)
	{
		solution[i] = sol_rev_skal[i + 2];
	}
	Norm2(solution, m);

	cout << endl;
	cout << "4) Противоположная граница спектра: " << sol_rev_skal[1] << endl;
	cout << "Соответствующий ему собственный вектор: " << endl;
	for (int i = 0; i < m; i++)
	{
		cout << setprecision(16) << solution[i] << " ";
	}
	cout << endl << endl;

	// Методом Виландта найдём собственное число и его вектор, используя степенной ход
	Vil(A, sol_vil, m);
	int k = sol_vil[0];
	cout << setprecision(16) << "5) Изолированное собственное число методом Виландта: " << sol_vil[k - 1] << endl << endl;

	//Уточнение Эйткена
	double sol_eyt = Eytken(sol_vil[k - 1], sol_vil[k - 2], sol_vil[k - 3]);
	cout << setprecision(16) << "6) Уточнение Эйткена изолированного собственного числа: " << sol_eyt << endl;


	return 0;
}

void Jacobi(double** Ak, double** X, int m)
{
	double** Ai = new double* [m]; // создаём массивы
	for (int i = 0; i < m; i++)
		Ai[i] = new double[m] {0};

	double** V = new double* [m];
	for (int i = 0; i < m; i++)
		V[i] = new double[m] {0};

	double** Xi = new double* [m]; // создаём массивы
	for (int i = 0; i < m; i++)
		Xi[i] = new double[m] {0};



	double maxaij = abs(Ak[0][1]);//A
	int maxi = 0;
	int maxj = 1;
	//double phi;
	double c, s;
	double d;

	//double cx, sx;
	//double dx;

	bool Ans = true;
	while (Ans == true)
	{
		for (int i = 0; i < m - 1; i++)
		{
			for (int j = i + 1; j < m; j++)
			{
				if (abs(Ak[i][j]) > maxaij)
				{
					maxaij = abs(Ak[i][j]);
					maxi = i;
					maxj = j;
				}
			}
		}
		//cout << "maxaij =" << maxaij << endl;
		if (abs(maxaij) < e)
		{
			Ans = false;
			break;
		}

		d = sqrt(pow(Ak[maxi][maxi] - Ak[maxj][maxj], 2) + 4 * pow(Ak[maxi][maxj], 2));
		c = sqrt(0.5 * (1 + abs(Ak[maxi][maxi] - Ak[maxj][maxj]) / d));
		int sgn;
		if ((Ak[maxi][maxj] * (Ak[maxi][maxi] - Ak[maxj][maxj])) == 0)
			sgn = 0;
		else if ((Ak[maxi][maxj] * (Ak[maxi][maxi] - Ak[maxj][maxj])) > 0)
			sgn = 1;
		else sgn = -1;
		s = sgn * sqrt(0.5 * (1 - abs(Ak[maxi][maxi] - Ak[maxj][maxj]) / d));

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++)
			{
				if ((i != maxi) && (i != maxj) && (i != maxj) && (j != maxi))
				{
					Ai[i][j] = Ak[i][j];

				}
				else if ((i != maxi) && (i != maxj))
				{
					Ai[i][maxi] = c * Ak[i][maxi] + s * Ak[i][maxj];
					Ai[maxi][i] = Ak[i][maxi];
					Ai[i][maxj] = (-1) * s * Ak[i][maxi] + c * Ak[i][maxj];
					Ai[maxj][i] = Ak[i][maxj];

				}
				
			}
		}
		Ai[maxi][maxi] = pow(c, 2) * Ak[maxi][maxi] + 2 * c * s * Ak[maxi][maxj] + pow(s, 2) * Ak[maxj][maxj];
		Ai[maxj][maxj] = pow(s, 2) * Ak[maxi][maxi] - 2 * c * s * Ak[maxi][maxj] + pow(c, 2) * Ak[maxj][maxj];
		Ai[maxi][maxj] = 0;
		Ai[maxj][maxi] = 0;

		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++)
			{
				if (i != maxi && i != maxj)
					V[i][i] = 1;
				if (i != maxi && i != maxj && j != maxi && j != maxj)
					V[i][j] = 0;
			}

		V[maxi][maxi] = c;
		V[maxj][maxj] = c;
		V[maxi][maxj] = (-1)*s;
		V[maxj][maxi] = s;

		MxM(X, V, Xi, m);


		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++)
			{
				Ak[i][j] = Ai[i][j];
			}
		}


		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++)
			{
				X[i][j] = Xi[i][j];
			}

		maxaij = abs(Ak[0][1]);
		maxi = 0;
		maxj = 1;

	}
}

void Norm2(double* X, int n)
{
	double summ=0;
	for (int i = 0; i < n; i++)
	{
		summ += pow(X[i],2);
	}
	summ = sqrt(summ);
	for (int i = 0; i < n; i++)
		X[i] = X[i] / summ;
}

void MxM(double** A, double** B, double** X, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				X[i][j] += A[i][k] * B[k][j];
}

void MxN(double** A, double* Y, double* Y_new, int n)
{
	for (int i = 0; i < n; i++)
		Y_new[i] = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Y_new[i] += A[i][j] * Y[j];
}

double Apos_est(double** A, double* Y, double lambda, int n)
{
	double* Y_new = new double[n];
	MxN(A, Y, Y_new, n);
	for (int i = 0; i < n; i++)
		Y[i] = lambda * Y[i];
	for (int i = 0; i < n; i++)
		Y_new[i] = Y_new[i] - Y[i];
	double YY_new = 0, YY = 0;
	for (int i = 0; i < n; i++)
	{
		YY_new += Y_new[i] * Y_new[i];
		YY += Y[i] * Y[i];
	}
	return sqrt(YY_new) / sqrt(YY);
}

void rev_Skalar(double** A, double* X, double lambda, int n)
{
	double rev_lambda = 0;
	double** B = new double* [n];
	for (int i = 0; i < n; i++)
		B[i] = new double[n] {0};

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (i == j)
				B[i][j] = A[i][j] - lambda;
			else
				B[i][j] = A[i][j];
		}
	Skalar(B, X, n);
	rev_lambda = X[1] + lambda;
	X[1] = rev_lambda;
}

void Skalar(double** A, double* sol, int n)
{
	int count = 0;
	double* Y_0 = new double[n];
	double* Y = new double[n] {0};
	double Apos = 1, lambda = 0;
	for (int i = 0; i < n; i++)
		Y_0[i] = 1;
	while (Apos > e)
	{
		MxN(A, Y_0, Y, n);
		double YY_0 = 0, Y0Y0 = 0;
		for (int i = 0; i < n; i++)
		{
			YY_0 += Y[i] * Y_0[i];
			Y0Y0 += Y_0[i] * Y_0[i];
		}
		lambda = YY_0 / Y0Y0;
		for (int i = 0; i < n; i++)
		{
			Y_0[i] = Y[i];
		}
		Apos = Apos_est(A, Y, lambda, n);
		count++;
	}
	sol[0] = count;
	sol[1] = lambda;
	for (int i = 2; i < n + 2; i++)
		sol[i] = Y[i - 2];
}

void Degree(double** A, double* sol, int n)
{
	int count = 0;
	double* Y_0 = new double[n];
	double* Y = new double[n];
	double Apos = 1, lambda = 0;
	for (int i = 0; i < n; i++)
		Y_0[i] = 1;

	while (Apos > e1)
	{
		MxN(A, Y_0, Y, n);
		lambda = Y[0] / Y_0[0];
		//cout << "lambda " << lambda << endl;
		for (int i = 0; i < n; i++)
			Y_0[i] = Y[i];
		Apos = Apos_est(A, Y, lambda, n);
		count++;
	}
	sol[0] = count;
	sol[1] = lambda;
	for (int i = 2; i < n + 2; i++)
		sol[i] = Y[i - 2];
}

void Vil(double** A, double* sol, int n)
{
	double lambda_0 = 1, lambda_k = 0, mu = 0, k = 1;
	int count = 0;

	double* Y_0 = new double[n];
	double* Y = new double[n];
	double** W = new double* [n];
	for (int i = 0; i < n; i++)
		W[i] = new double[n+1] {0};

	for (int i = 0; i < n; i++)
		Y_0[i] = 1;
	while (k > e1)
	{
		for (int i = 0; i < n; i++)
		{
			Y_0[i] = 1;
			for (int j = 0; j < n; j++)
			{
				if (i == j)
					W[i][j] = A[i][j] - lambda_0;
				else
					W[i][j] = A[i][j];
			}
		}
		for (int i = 0; i < n; i++)
		{
			W[i][n] = Y_0[i];
		}
		Gauss(W, Y, n, n+1);

		double YY0 = 0, Y0Y0=0;
		for (int i = 0; i < n; i++)
		{
			YY0 += Y[i] * Y_0[i];
			Y0Y0 += Y_0[i] * Y_0[i];
		}

		mu = YY0 / Y0Y0;
		lambda_k = 1 / (mu)+lambda_0;
		k = abs(lambda_k - lambda_0);
		for (int i = 0; i < n; i++)
		{
			Y_0[i] = Y[i];
		}
		lambda_0 = lambda_k;
		sol[count+1]=lambda_k;
		count++;
	}
	sol[0] = count;
}

void Gauss(double** A, double* X, int m, int n) // с выбором главного элемента по столбцу и строке
{
	double** A1 = new double* [m];
	for (int i = 0; i < m; i++)
		A1[i] = new double[n];

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			A1[i][j] = A[i][j];
	}

	double* Xi = new double[m];
	int* iX = new int[m];
	for (int i = 0; i < m; i++)
	{
		iX[i] = i;
	}
	int t1 = 0, t2 = 0, index;
	double Akt = -1;
	for (int k = 0; k < m - 1; k++)
	{
		t1 = k, t2 = k;
		for (int i = k; i < m; i++)
		{
			for (int j = k; j < m; j++)
			{
				if (abs(A1[i][j]) > Akt)
				{
					Akt = A1[i][j];
					t1 = i;
					t2 = j;
				}
			}
		}

		if (t2 != k)
		{
			index = iX[t2];
			iX[t2] = iX[k];
			iX[k] = index;
		}

		for (int i = 0; i < m; i++)
		{
			Akt = A1[i][k];
			A1[i][k] = A1[i][t2];
			A1[i][t2] = Akt;
		}

		for (int i = 0; i < n; i++)
		{
			Akt = A1[k][i];
			A1[k][i] = A1[t1][i];
			A1[t1][i] = Akt;
		}
		for (int j = n - 1; j >= k; j--)
		{
			if (A1[k][k] != 0)
			{
				A1[k][j] = A1[k][j] / A1[k][k];
			}
		}

		for (int i = k + 1; i < m; i++)
		{
			for (int j = n - 1; j >= k; j--)
			{
				A1[i][j] = A1[i][j] - A1[k][j] * A1[i][k];
			}
		}
		Akt = -1;
	}

	for (int j = n - 1; j >= 0; j--)
	{
		A1[m - 1][j] = A1[m - 1][j] / A1[m - 1][m - 1];
	}

	double temp = 0;
	for (int i = 0; i < m; i++)
	{
		X[i] = 0;
		Xi[i] = 0;
	}

	for (int k = m - 1; k >= 0; k--)
	{
		for (int j = k; j < n; j++)
		{
			temp = temp + Xi[j] * A1[k][j];
		}
		Xi[k] = A1[k][n - 1] - temp;
		temp = 0;
	}
	for (int i = 0; i < m; i++)
	{
		X[iX[i]] = Xi[i];
	}
}

double Eytken(double a, double b, double c)
{
	return (a * c - pow(b, 2)) / (a - 2 * b + c);
}