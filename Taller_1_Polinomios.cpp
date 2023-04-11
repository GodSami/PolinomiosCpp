#include <iostream>
#include <complex>
#include <vector>
#define M_PI 3.1415926535
std::vector<std::complex<double>> durand_kerner(std::vector<std::complex<double>> coefficients) {
    int degree = coefficients.size() - 1;
    std::vector<std::complex<double>> roots(degree);
    for (int i = 0; i < degree; i++) {
        roots[i] = std::polar(1.0, 2 * M_PI * i / degree);
    }
    const int max_iterations = 1000;
    const double tolerance = 1e-6;
    for (int iteration = 0; iteration < max_iterations; iteration++) {
        bool converged = true;
        std::vector<std::complex<double>> new_roots(degree);
        for (int i = 0; i < degree; i++) {
            std::complex<double> numerator = coefficients[degree];
            std::complex<double> denominator = 1.0;
            for (int j = 0; j < degree; j++) {
                numerator = numerator * roots[i] + coefficients[degree - j - 1];
                if (j != i) {
                    denominator *= roots[i] - roots[j];
                }
            }
            new_roots[i] = roots[i] - numerator / denominator;
            if (std::abs(new_roots[i] - roots[i]) > tolerance) {
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
    std::cout << "Ingresa el grado del polinomio: ";
    int degree;
    std::cin >> degree;
    std::vector<std::complex<double>> coefficients(degree + 1);
    for (int i = 0; i <= degree; i++) {
        std::cout << "Ingresa el coeficiente de x^" << degree - i << ": ";
        std::cin >> coefficients[i];
    }
    std::vector<std::complex<double>> roots = durand_kerner(coefficients);
    for (const auto& root : roots) {
        std::cout << "La raiz es: " << root << std::endl;
    }
    return 0;
}
