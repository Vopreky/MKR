#include <iostream>
using namespace std;

//граничные условия
//(x,y) из границы([0, 1] X [0,1])
float g(float x, float y){
	return x*y;
}
//функиця f
//(x,y) из ([0, 1] X [0,1])
float f(float x, float y){
	return 0;
}
void print(float * A, float * b, int size){
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
			cout << A[i * size + j] << '\t';
		}
		cout << '\t' << b[i] << endl;
	}
}
void solve(float * A, float * b, float * x, int size){
	int * tmp = new int[size];
	for (int i = 0; i < size; i++){
		tmp[i] = i;
	}
	for (int i = 0; i < size; i++){
		int ind = i;
		for (int j = i + 1; j < size; j++){
			if (abs(A[i * size + tmp[ind]]) < abs(A[i * size + tmp[j]])){
				ind = j;
			}
		}
		if (ind != i){
			int t = tmp[i];
			tmp[i] = tmp[ind];
			tmp[ind] = t;
		}
		for (int j = 0; j < size; j++){
			if (j == i)
				continue;
			float coef = A[j * size + tmp[i]] / A[i * size + tmp[i]];
			b[j] -= b[i] * coef;
			for (int k = 0; k < size; k++){
				A[j * size + k] -= A[i * size + k] * coef;
			}
		}
	}
	for (int i = 0; i < size; i++){
		x[tmp[i]] = b[i] / A[i * size + tmp[i]];
	}
	delete[] tmp;
}
int main(int argc, char ** argv){
	const float d11 = 1, d12 = 0, d22 = 1;
	//размеры сетки
	const int sizeX = 17;
	const int sizeY = 17;
	//матрица системы
	float A[(sizeX * sizeY) * (sizeX * sizeY)] = {};
	//правая часть
	float b[sizeX * sizeY] = {};
	//тут будет ответ
	float x[sizeX * sizeY] = {};

	const float hx = 1.0 / (sizeX - 1);
	const float hy = 1.0 / (sizeY - 1);

	for (int i = 0; i < sizeX; i++){
		for (int j = 0; j < sizeY; j++){
			if (i == 0){
				b[i * sizeY + j] = g(hx * i, hy * j);
				A[(i * sizeY + j) * (sizeX * sizeY) + (i * sizeY + j)] = 1;
			}
			else if (i == sizeX - 1){
				b[i * sizeY + j] = g(hx * i, hy * j);
				A[(i * sizeY + j) * (sizeX * sizeY) + (i * sizeY + j)] = 1;
			}
			else if (j == 0){
				b[i * sizeY + j] = g(hx * i, hy * j);
				A[(i * sizeY + j) * (sizeX * sizeY) + (i * sizeY + j)] = 1;
			}
			else if (j == sizeY - 1){
				b[i * sizeY + j] = g(hx * i, hy * j);
				A[(i * sizeY + j) * (sizeX * sizeY) + (i * sizeY + j)] = 1;
			}
			else{
				A[(i * sizeY + j) * (sizeX * sizeY) + (i * sizeY + j)] = d11 * -2.0 / hx / hx + d22 * -2.0 / hy / hy;
				A[(i * sizeY + j) * (sizeX * sizeY) + (i * sizeY + j - 1)] = d22 * 1.0 / hy / hy;
				A[(i * sizeY + j) * (sizeX * sizeY) + (i * sizeY + j + 1)] = d22 * 1.0 / hy / hy;
				A[(i * sizeY + j) * (sizeX * sizeY) + ((i - 1) * sizeY + j)] = d11 * 1.0 / hx / hx;
				A[(i * sizeY + j) * (sizeX * sizeY) + ((i + 1) * sizeY + j)] = d11 * 1.0 / hx / hx;
				A[(i * sizeY + j) * (sizeX * sizeY) + ((i - 1) * sizeY + j - 1)] = d12 * 1.0 / hx / hy;
				A[(i * sizeY + j) * (sizeX * sizeY) + ((i - 1) * sizeY + j + 1)] = d12 * -1.0 / hx / hy;
				A[(i * sizeY + j) * (sizeX * sizeY) + ((i + 1) * sizeY + j - 1)] = d12 * -1.0 / hx / hy;
				A[(i * sizeY + j) * (sizeX * sizeY) + ((i + 1) * sizeY + j + 1)] = d12 * 1.0 / hx / hy;
				b[i * sizeY + j] = f(hx * i, hy * j);
			}
		}
	}
	print(A, b, sizeX * sizeY);
	cout << endl << "ответ: " << endl;
	solve(A, b, x, sizeX * sizeY);
	for (int i = 0; i < sizeX * sizeY; i++){
		cout << "u[" << i / sizeY << ", " << i % sizeY << "] = u(" << i / sizeY * hx << ", " << i % sizeY * hy << ") = " << x[i] << endl;
	}
}
