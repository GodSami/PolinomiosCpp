#include <iostream>
#include <complex>
#include <vector>
#define M_PI 3.1415926535

using namespace std;

vector<complex<double>> durand_kerner(vector<complex<double>> coefficients) {
    int degree = coefficients.size() - 1;
    vector<complex<double>> roots(degree);
    for (int i = 0; i < degree; i++) {
        roots[i] = polar(1.0, 2 * M_PI * i / degree);
    }
    const int max_iterations = 1000;
    const double tolerance = 1e-6;
    for (int iteration = 0; iteration < max_iterations; iteration++) {
        bool converged = true;
        vector<complex<double>> new_roots(degree);
        for (int i = 0; i < degree; i++) {
            complex<double> numerator = coefficients[degree];
            complex<double> denominator = 1.0;
            for (int j = 0; j < degree; j++) {
                numerator = numerator * roots[i] + coefficients[degree - j - 1];
                if (j != i) {
                    denominator *= roots[i] - roots[j];
                }
            }
            new_roots[i] = roots[i] - numerator / denominator;
            if (abs(new_roots[i] - roots[i]) > tolerance) {
                converged = false;
            }
        }
        roots = new_roots;
        if (converged) {
            break;
        }
    }
    return roots;
}

int main() {
    cout << "Ingresa el grado del polinomio: ";
    int degree;
    cin >> degree;
    vector<complex<double>> coefficients(degree + 1);
    for (int i = 0; i <= degree; i++) {
        cout << "Ingresa el coeficiente de x^" << degree - i << ": ";
        cin >> coefficients[i];
    }
    vector<complex<double>> roots = durand_kerner(coefficients);
    for (const auto& root : roots) {
        cout << "La raiz es: " << root << endl;
    }
    return 0;
}
