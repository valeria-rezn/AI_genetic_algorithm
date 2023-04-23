#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <ctime>   
#include <cstdlib>
#include <bitset>
#include <algorithm>
using namespace std;
//test
double minim(vector <double>& fit)
{
	double min = fit[0];
	for (int i = 1; i < fit.size(); i++)
	{
		if (min > fit[i])
			min = fit[i];
	}
	return min;
}
void  maxim(vector <double>& fit, vector <double>& x, vector <double>& y, vector <double>& x_m, vector <double>& y_m, vector <double>& fit_m)
{
	double max = fit[0];
	int j = 0;
	for (int i = 1; i < fit.size(); i++)
	{
		if (max < fit[i])
		{
			max = fit[i];
			j = i;
		}
	}
	x_m.push_back(x[j]);
	y_m.push_back(y[j]);
	fit_m.push_back(max);

}
void fitness_function(vector <double>&x, vector <double>&y, vector <double>&f_ind, vector <double>&f_avg, int N)
{
	double fit_sum = 0;
	for (int i = 0; i < N; i++)
	{
		f_ind.push_back(2*x[i] + y[i]);
		fit_sum += f_ind[i];
	}
	f_avg.push_back(fit_sum / N);//average fitness value
}
//normalization fitness function for roulette
void normalization(int N, vector <double> &fit_i,double min)
{
	double f_sum = 0;
	for (int i = 0; i < N; i++)
	{
		fit_i[i] += 2 * abs(min);
		f_sum += fit_i[i];
	}
	for (int i = 0; i < N; i++)
	{
		fit_i[i] = fit_i[i]/f_sum;
	}
}
//choose random ind from population
void roulette(int N, vector <double> &f_ind, vector<int>& id)
{
	for(int j=0;j<N;j++)
	{
		double r =0.99*(double)rand() / (double)RAND_MAX+0.01;
		int i = 0;
		while (r > 0 && i<N)
		{
			r -= f_ind[i];
			i++;
		}
		cout << endl;
		id.push_back(i-1);
		//cout <<j <<"    id = " << id[j] << endl;
	}
}
void chooseParents(vector <double> &x, vector <double> &x1, vector <double>&y, vector <double>&y1,vector<int>&id, vector<double>&f_ind,int N)
{
	for (int i = 0; i < N; i++)
	{
		x1.push_back(x[id[i]]);
		//cout <<i <<"    x1 = " << x1[i]<<endl;
		y1.push_back(y[id[i]]);
		//cout <<i<< "    y1 = " << y1[i] << endl;
	}
	//cout << "DEATH" << endl;
	x = x1;
	y = y1;
}
void codingNumbers(vector <double>&x,vector <int>&x_10, int a, double h,int N)
{
	for (int i = 0; i < N; i++)
	{
		x_10.push_back(roundf((x[i] - a) / h));
	}

}

void toBinary(vector <int> &x_10,vector <string> &x_bi, int l,int N )
{
	for (int i = 0; i < N; i++)
	{
		string r_t;
		while (x_10[i] != 0)
		{
			r_t = (x_10[i] % 2 == 0 ? "0" : "1") + r_t;
			x_10[i] /= 2;
		}
		string result;
		int dif = l - r_t.length();
		while (dif != 0)
		{
			result += "0";
			dif--;
		}
		x_bi.push_back(result+r_t);
	}
}
void crossover(int N,int l,vector <string> &x_bi)
{
	double properity = 0.95;
	int k = 0;
	vector <string> temp;
	temp.push_back(x_bi[0]);
	for (int i = 0; i < N; i += 2)
	{
		double lambda = rand() / RAND_MAX;
		if (lambda < properity)
		{
			k = 1+rand() % (l - 2);
			temp[0] = x_bi[i];
			for (int j = k; j < l; j++)
			{
				x_bi[i][j] = x_bi[i + 1][j];
				x_bi[i + 1][j] = temp[0][j];
			}
		}
	}
}
void mutation(int N, int l, vector <string>&x_bi)
{
	
	double properity = 0.1;
	int k = 0;
	for (int i = 0; i < N; i += 2)
	{
		double lambda = rand() / RAND_MAX;
		if (lambda < properity)
		{
			k = rand() % (l - 1);
			x_bi[i][k] = (x_bi[i][k] == '0' ? '1' : '0');
		}

	}

}
void fromBinaryToCoded(vector <string>&x_bi, vector <int>&x_10,int N,int l)
{
	vector<int> temp(l, 0);
	int base = 1;
	for (int i = 0; i < N; i++)
	{
		int sum = 0;
		base = 1;
		int power = l - 1;
		for (int j = 0; j <l; j++)
		{
				if (x_bi[i][j] == '1')
				{
					sum +=pow(2,power);
				}
				power--;
		}
		x_10[i] = sum;
		temp.resize(l,0);
	}
}
void fromCodedToStart(vector <int>& x_10, vector <double>& x, int N,int a,double h)
{
	
	for (int i = 0; i < N; i++)
	{
		x[i] = a + x_10[i] * h;
	}
	
}
void showPopulation(int N, int num, vector <double>& x, vector < double > & y)
{
	for(int i=0;i<N;i++)
	cout << num<< "   X = " << x[i] << "   Y = " << y[i] << endl;
	cout << endl << endl;
}
int main() {
	//оголошення вхідних даних
	double a, b, c, d;
	int q, N, itteration;
	cout << "enter value of LEFT border for X: ";
	cin >> a;
	cout << endl << "enter value of RIGHT border for X: ";
	cin >> b;
	cout << endl << "enter value of LEFT border for Y: ";
	cin >> c;
	cout << endl << "enter value of RIGHT border for Y: ";
	cin >> d;
	cout << endl << "enter number of decimal places: ";
	cin >> q;
	cout << endl << "enter population size: ";
	cin >> N;
	cout << endl << "enter number of itteration";
	cin >> itteration;
	//розрахунок необхідних вхідних даних
	int length_x = roundf(log2((b - a) * pow(10, q) + 1));
	int length_y = roundf(log2((d - c) * pow(10, q) + 1));
	double N_x = pow(2, length_x) - 1;
	double N_y = pow(2, length_y) - 1;
	double h_x = (b - a) / N_x;
	double h_y = (d - c) / N_y;
	int num_pop = 0;//number of population
	double minimal = 0;
	vector <double> X;
	vector <double> Y;
	vector <double> fit_ind;//individual fitness
	vector <double> fitness_avarage;
	vector <double> X_1;
	vector <double> Y_1;
	vector <int> X_10;
	vector <int> Y_10;
	vector <string> X_bi;
	vector <string> Y_bi;
	vector <int> id;
	vector < double> max_fit;
	vector <double> max_x;
	vector <double> max_y;
	srand(time(NULL));
	for (int i = 0; i < N; i++)
	{
		X.push_back(a + (rand() / (RAND_MAX / (b - a))));
		Y.push_back(c + (rand() / (RAND_MAX / (d - c))));
	}
	showPopulation(N, num_pop, X, Y);
	//start of terriful things
	for (int j = 0; j < itteration; j++)
	{
		fitness_function(X, Y, fit_ind, fitness_avarage, N);
		//cout << "f_norm" << endl;
		minimal = minim(fit_ind);
		maxim(fit_ind, X, Y, max_x, max_y, max_fit);
		//cout << "max" << endl;
		normalization(N, fit_ind, minimal);
		//cout << "norm norm"<<endl;
		roulette(N, fit_ind, id);
		//cout << "roul ok" << endl;
		chooseParents(X, X_1, Y, Y_1, id, fit_ind, N);
		//cout << "par norm" << endl;
		codingNumbers(X, X_10, a, h_x, N);
		//cout << "cod_x" << endl;
		codingNumbers(Y, Y_10, c, h_y, N);
		//cout << "cod_y" << endl;
		toBinary(X_10, X_bi, length_x, N);
		//cout << "bi x" << endl;
		toBinary(Y_10, Y_bi, length_y, N);
		//cout << "bi y" << endl;
		crossover(N, length_x, X_bi);
		//cout << "cr x" << endl;
		crossover(N, length_y, Y_bi);
		//cout << "cr y" << endl;
		mutation(N, length_x, X_bi);
		//cout << "m x" << endl;
		mutation(N, length_y, Y_bi);
		//cout << "m y" << endl;
		fromBinaryToCoded(X_bi, X_10, N, length_x);
		//cout << "fbtc x" << endl;
		fromBinaryToCoded(Y_bi, Y_10, N, length_y);
		//cout << "fbtc y" << endl;
		fromCodedToStart(X_10, X, N, a, h_x);
		//cout << "fcts x" << endl;
		fromCodedToStart(Y_10, Y, N, c, h_y);
		//cout << "fcts y" << endl;
		num_pop++;
		//cout << "NUMBER OF POP: " << num_pop << endl << endl;
		showPopulation(N, num_pop, X, Y);
		//clear all temp
		X_1.clear();
		Y_1.clear();
		X_10.clear();
		Y_10.clear();
		X_bi.clear();
		Y_bi.clear();
		fit_ind.clear();
		id.clear();
		cout << endl << endl << endl << endl;
	}
	double max = 0;
	int j = 0;
	for (int i = 0; i < max_fit.size(); i++)
	{
		cout << max_fit[i] << endl;
		if (max < max_fit[i])
		{
			max = max_fit[i];
			j = i;
		}
	}
	cout<<"max value of fitness function is "<<max<<" in dote with X = "<<max_x[j]<<" and Y = "<<max_y[j];
}

