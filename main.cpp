#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double g(double x, double u) {
    return 1000 * (1 - u * u) * u - x;
}

vector<double> Runge_Kutta4(double a, double b, double x0, double u0, int N) {
    double h = (b - a) / N;
    vector<double> t;
    t.resize(N);
    t[0] = a;
    for (int i = 1; i < N; i++)
        t[i] = t[i-1] + h;
    vector<double> x;
    x.resize(N);
    vector<double> u;
    u.resize(N);
    x[0] = x0;
    u[0] = u0;
    double k1, k2, k3, k4, l1, l2, l3, l4;
    for (int i = 1; i < N; i++) {
        k1 = h * u[i - 1];
        l1 = h * g( x[i - 1], u[i - 1]);
        k2 = h * (u[i - 1] + l1 / 2);
        l2 = h * g(x[i - 1] + k1 / 2, u[i - 1] + l1 / 2);
        k3 = h * (u[i - 1] + l2 / 2);
        l3 = h * g( x[i - 1] + k2 / 2, u[i - 1] + l2 / 2);
        k4 = h * (u[i - 1] + l3);
        l4 = h * g( x[i - 1] + k3 / 2, u[i - 1] + l3 / 2);
        x[i] = x[i - 1] + 1. / 6 * (k1 + 2 * (k2 + k3) + k4);
        u[i] = u[i - 1] + 1. / 6 * (l1 + 2 * (l2 + l3) + l4);

    }
    cout<<t[1]<<endl;
    return x;

}
vector<double> RK (double a, double b, double x0, double u0, int N){
    double h = (b-a)/N;
    vector<double> t;
    t.resize(N);
    t[0] = a;
    for (int i = 1; i < N; i++)
        t[i] += i * h;
    vector<double> x;
    x.resize(N);
    vector<double> u;
    u.resize(N);
    x[0] = x0;
    u[0] = u0;
    for (int i = 1; i < N; i++) {
        x[i] = x[i-1] + h*u[i-1] + (h*h*u[i-1]/2);
        u[i] = u[i-1] + h*g(x[i-1],u[i-1]);

    }
    return x;
}

int main() {
    int N = 1000000;
    double a = 0.;
    double b = 1000.;
    double x0 = 0.;
    double u0 = 0.001;
    vector<double> res = Runge_Kutta4(a, b, x0, u0, N);
    //vector<double> res1 = RK(a, b, x0, u0, N);
    for(int i = 1;i < 100;i++)
        cout<<res[i*10000]<<endl;

}
