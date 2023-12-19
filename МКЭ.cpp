#include <unistd.h>
#include <iostream>
#include <list>
#include <cmath>
//#include <SDL2/SDL.h>
const float d11 = 1, d12 = 0, d21 = 0, d22 = 1;

using namespace std;
int width, height;


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

float g(float x, float y){
	return cos(x)*sin(y);
}
float g(float * x){
	return g(x[0], x[1]);
}
float f(float x, float y){
	return -2*cos(x)*sin(y);
}
float f(float * x){
	return f(x[0], x[1]);
}

class FEM{
	int sizeX, sizeY;
	float * pointArray;
	float * values;
	list<int> * con;
public:
	FEM(int sizeX, int sizeY){
		this->sizeX = sizeX;
		this->sizeY = sizeY;
		int N = sizeX * sizeY;
		values = new float[N];
		pointArray = new float[sizeX * sizeY * 2];
		con = new list<int>[sizeX * sizeY];
		for (int i = 0; i < sizeX; i++)
		for (int j = 0; j < sizeY; j++){
			pointArray[2 * (j * sizeX + i) + 0] = i / (float)(sizeX - 1);
			pointArray[2 * (j * sizeX + i) + 1] = j / (float)(sizeY - 1);
			int ci, cj;
			
			ci = i - 1;
			cj = j - 1;
			if (ci >= 0 && ci < sizeX && cj >= 0 && cj < sizeY){
				con[(j * sizeX + i)].push_back(cj * sizeX + ci);
			}
			ci = i;
			cj = j - 1;
			if (ci >= 0 && ci < sizeX && cj >= 0 && cj < sizeY){
				con[(j * sizeX + i)].push_back(cj * sizeX + ci);
			}
			ci = i - 1;
			cj = j;
			if (ci >= 0 && ci < sizeX && cj >= 0 && cj < sizeY){
				con[(j * sizeX + i)].push_back(cj * sizeX + ci);
			}
			ci = i + 1;
			cj = j;
			if (ci >= 0 && ci < sizeX && cj >= 0 && cj < sizeY){
				con[(j * sizeX + i)].push_back(cj * sizeX + ci);
			}
			ci = i;
			cj = j + 1;
			if (ci >= 0 && ci < sizeX && cj >= 0 && cj < sizeY){
				con[(j * sizeX + i)].push_back(cj * sizeX + ci);
			}
			ci = i + 1;
			cj = j + 1;
			if (ci >= 0 && ci < sizeX && cj >= 0 && cj < sizeY){
				con[(j * sizeX + i)].push_back(cj * sizeX + ci);
			}
		}
		for (int i = 0; i < N; i++){
			for (auto cj = con[i].begin(); cj != con[i].end(); cj++){
				int j = *cj;
				cout << i << " -> " << j << endl;
			}
		}
	}
	/*
	void draw(SDL_Renderer * rend) const{
		SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
		for (int i = 0; i < sizeY*sizeX; i++){
			float x, y;
			x = pointArray[2*(i) + 0];
			y = pointArray[2*(i) + 1];
			int px, py;
			px = x * (width - 1);
			py = y * (height - 1);
			SDL_RenderDrawLine(rend, px - 21, py, px + 21, py);
			SDL_RenderDrawLine(rend, px, py - 21, px, py + 21);
			for (auto j = con[i].begin(); j != con[i].end(); j++){
				float cx, cy;
				cx = pointArray[2*(*j) + 0];
				cy = pointArray[2*(*j) + 1];
				int pcx, pcy;
				pcx = cx * (width - 1);
				pcy = cy * (height - 1);
				SDL_RenderDrawLine(rend, px, py, pcx, pcy);

			}
		}
		return;
	}
	*/
	void findNN(int ind1, int ind2, int * na, int * count){
		list<int> & l1 = con[ind1];
		list<int> & l2 = con[ind2];
		int cn = 0;
		for (auto i = l1.begin(); i != l1.end(); i++){
			for (auto j = l2.begin(); j != l2.end(); j++){
				if (*i == *j){
					na[cn++] = *i;
				}
			}
		}
		*count = cn;
	}
	void loadPoint(int ind, float * vec) const{
		for (int i = 0; i < 2; i++){
			vec[i] = pointArray[ind * 2 + i];
		}
	}
	float bf(const float * x, int i1, int i2, int i3) const{
		float p1[2] = {};
		float p2[2] = {};
		float p3[2] = {};
		loadPoint(i1, p1);
		loadPoint(i2, p2);
		loadPoint(i3, p3);
		float A[9] = {1, 0, 0, 1, 0, 0, 1, 0, 0};
		float c[3] = {};
		float b[3] = {1, 0, 0};
		for (int j = 0; j < 2; j++){
			A[0 * 3 + 1 + j] = p1[j];
		}
		for (int j = 0; j < 2; j++){
			A[1 * 3 + 1 + j] = p2[j];
		}
		for (int j = 0; j < 2; j++){
			A[2 * 3 + 1 + j] = p3[j];
		}
		solve(A, b, c, 3);
		float c1 = c[0], c2 = c[1], c3 = c[2];
		/*
		float p1[2] = {};
		float p2[2] = {};
		float p3[2] = {};
		loadPoint(i1, p1);
		loadPoint(i2, p2);
		loadPoint(i3, p3);
		cout << p1[0] << ", " << p1[1] << endl;
		cout << p2[0] << ", " << p2[1] << endl;
		cout << p3[0] << ", " << p3[1] << endl;
		if (p3[0] != p2[0]){
			c3 = ((p3[1] - p1[1])*(p2[0] - p1[0]) - (p2[1] - p1[1])*(p3[0] - p1[0])) / (p3[0] - p2[0]);
		}
		else{
			c3 = 0;
		}
		if (p3[1] != p2[1]){
			c2 = ((p3[0] - p1[0])*(p2[1] - p1[1]) - (p2[0] - p1[0])*(p3[1] - p1[1])) / (p3[1] - p2[1]);
		}
		else{
			c2 = 0;
		}
		c1 = 1 - p1[0] * c2 - p1[1] * c3;
		*/
		cout << c1 << " " << c2 << " " << c3 << endl;
		return (c1 + c2 * x[0] + c3 * x[1]);

	}
	float bd(int dim, int i1, int i2, int i3) const{
		float p1[2] = {};
		float p2[2] = {};
		float p3[2] = {};
		loadPoint(i1, p1);
		loadPoint(i2, p2);
		loadPoint(i3, p3);
		float A[9] = {1, 0, 0, 1, 0, 0, 1, 0, 0};
		float c[3] = {};
		float b[3] = {1, 0, 0};
		for (int j = 0; j < 2; j++){
			A[0 * 3 + 1 + j] = p1[j];
		}
		for (int j = 0; j < 2; j++){
			A[1 * 3 + 1 + j] = p2[j];
		}
		for (int j = 0; j < 2; j++){
			A[2 * 3 + 1 + j] = p3[j];
		}
		solve(A, b, c, 3);
		float c1 = c[0], c2 = c[1], c3 = c[2];
		/*
		float c1 = 0, c2 = 0, c3 = 0;
		float p1[2] = {};
		float p2[2] = {};
		float p3[2] = {};
		loadPoint(i1, p1);
		loadPoint(i2, p2);
		loadPoint(i3, p3);
		if (p3[0] != p2[0]){
			c3 = ((p3[1] - p1[1])*(p2[0] - p1[0]) - (p2[1] - p1[1])*(p3[0] - p1[0])) / (p3[0] - p2[0]);
		}
		else{
			c3 = 0;
		}
		if (p3[1] != p2[1]){
			c2 = ((p3[0] - p1[0])*(p2[1] - p1[1]) - (p2[0] - p1[0])*(p3[1] - p1[1])) / (p3[1] - p2[1]);
		}
		else{
			c2 = 0;
		}
		c1 = 1 - p1[0] * c2 - p1[1] * c3;
		*/
		if (dim == 0)
			return c2;
		return c3;

	}
	float area (int i1, int i2, int i3){
		float * p1 = pointArray + i1 * 2;
		float * p2 = pointArray + i2 * 2;
		float * p3 = pointArray + i3 * 2;
		return area(p1, p2, p3);
	}
	float area(float * p1, float * p2, float * p3){
		float a, b, c, p;
		a = sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p2[1])*(p1[1] - p2[1]));
		b = sqrt((p2[0] - p3[0])*(p2[0] - p3[0]) + (p2[1] - p3[1])*(p2[1] - p3[1]));
		c = sqrt((p3[0] - p1[0])*(p3[0] - p1[0]) + (p3[1] - p1[1])*(p3[1] - p1[1]));
		p = (a + b + c) / 2;
		return sqrt(p * (p - a) * (p - b) * (p - c));
	}
	void setValues(){
		int N = sizeX * sizeY;
		float * A = new float[N * N];
		float * b = new float[N];
		for (int i = 0; i < N * N; i++){
			A[i] = 0;
		}
		for (int i = 0; i < N; i++){
			b[i] = 0;
		}
		for (int i = 0; i < N; i++){
			for (auto cj = con[i].begin(); cj != con[i].end(); cj++){
				int j = *cj;
				int na[2];
				int count;
				findNN(i, j, na, &count);
				/*
				cout << "nn(" << i << ", " << j << "):";
				for (int k = 0; k < count; k++){
					cout << " " << na[k];
				}
				cout << endl;
				*/
				float p1[2], p2[2], p3[2];
				loadPoint(i, p1);
				loadPoint(j, p2);
				for (int ck = 0; ck < count; ck++){
					int k = na[ck];
					loadPoint(k, p3);
					float cen[2];
					cen[0] = p1[0] / 3 + p2[0] / 3 + p3[0] / 3;
					cen[1] = p1[1] / 3 + p2[1] / 3 + p3[1] / 3;
					A[i * N + j] += area(p1, p2, p3) * d11 * (bd(0, i, j, k) * bd(0, j, i, k));
					A[i * N + j] += area(p1, p2, p3) * d12 * (bd(0, i, j, k) * bd(1, j, i, k));
					A[i * N + j] += area(p1, p2, p3) * d21 * (bd(1, i, j, k) * bd(0, j, i, k));
					A[i * N + j] += area(p1, p2, p3) * d22 * (bd(1, i, j, k) * bd(1, j, i, k));
					//b[i] += f(cen) * bf(cen, i, j, k) * area(p1, p2, p3) / 2;
					cout << "(" << i << ", " << j << ", " << k << ")" << endl;
					
					cout << "(" << i << ", " << i << ", " << j << " - " << k << ") / 2" << endl;
					A[i * N + i] += area(p1, p2, p3) * d11 * (bd(0, i, j, k) * bd(0, i, j, k)) / 2;
					A[i * N + i] += area(p1, p2, p3) * d12 * (bd(0, i, j, k) * bd(1, i, j, k)) / 2;
					A[i * N + i] += area(p1, p2, p3) * d21 * (bd(1, i, j, k) * bd(0, i, j, k)) / 2;
					A[i * N + i] += area(p1, p2, p3) * d22 * (bd(1, i, j, k) * bd(1, i, j, k)) / 2;
					cout << b[i] << " -> ";
					b[i] -= f(cen) * 0.333333333333 * area(p1, p2, p3) / 2; // (-16.0)
					cout << b[i] << endl;
					cout << "f(" << cen[0] << ", " << cen[1] << ") = " << f(cen) << endl;
					//b[i] += f(cen) * bf(cen, i, j, k) * area(p1, p2, p3) / 2;

					//cout << bf(cen, i, j, k) << endl;
				}

			}
			
		}
		for (int ix = 0; ix < sizeX; ix++){
			int ind1 = 0 * sizeX + ix;
			int ind2 = (sizeY - 1) * sizeX + ix;

			b[ind1] = g(ix / (float)(sizeX - 1), 0);
			b[ind2] = g(ix / (float)(sizeX - 1), 1);
			for (int j = 0; j < N; j++){
				A[(ind1) * N + j] = 0;
				A[(ind2) * N + j] = 0;
			}
			A[(ind1) * (N + 1)] = 1;
			A[(ind2) * (N + 1)] = 1;
		}
		cout << endl;
		for (int iy = 0; iy < sizeY; iy++){
			int ind1 = iy * sizeX + 0;
			int ind2 = iy * sizeX + sizeX - 1;

			b[ind1] = g(0, iy / (float)(sizeY - 1));
			b[ind2] = g(1, iy / (float)(sizeY - 1));
			for (int j = 0; j < N; j++){
				A[(ind1) * N + j] = 0;
				A[(ind2) * N + j] = 0;
			}
			A[(ind1) * (N + 1)] = 1;
			A[(ind2) * (N + 1)] = 1;
		}
		cout << endl;
		print(A, b, N);
		solve(A, b, values, N);
		for (int i = 0; i < N; i++){
			float x, y;
			int ix, iy;
			ix = i % sizeX;
			iy = i / sizeY;
			x = ix / (float)(sizeX - 1);
			y = iy / (float)(sizeY - 1);
			cout << i << ": u(" << x << ", " << y << ") = " << values[i] << endl;
		}
	}
	float errorNorm(){
		float err = 0;
		int N = sizeX * sizeY;
		for (int i = 0; i < N; i++){
			float x, y;
			int ix, iy;
			ix = i % sizeX;
			iy = i / sizeY;
			x = ix / (float)(sizeX - 1);
			y = iy / (float)(sizeY - 1);
			err = max(err, abs(values[i] - g(x, y)));
		}
		return err;
	}
};

int main(){
	FEM fem(25,25);
	float x1[2] = {0.5, 0.5};
	cout << endl;
	cout << endl;
	fem.setValues();
	cout << "error: " << fem.errorNorm() << endl;
	/*
	SDL_Init(SDL_INIT_VIDEO);
	SDL_DisplayMode DM;
	SDL_GetCurrentDisplayMode(0, &DM);
	width = DM.w;
	height = DM.h;
	SDL_Window * win = SDL_CreateWindow("windowTitle", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_WINDOW_FULLSCREEN);
	SDL_Renderer * rend = SDL_CreateRenderer(win, -1, 0);
	SDL_Event event;
	while (1){
		SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
		SDL_RenderClear(rend);
		fem.draw(rend);
		SDL_RenderPresent(rend);
		usleep(100000);
		while(SDL_PollEvent(&event)){
			if (event.type == SDL_QUIT){
				SDL_Quit();
				return 0;
			}
		}
	}
	*/

}
