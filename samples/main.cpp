// #include <glad/glad.h>
// #include <GLFW/glfw3.h>

#define DIFFEQ_FLOAT_PRECISION
#include "diffeq.h"

#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace DES;


int main() {

	auto t1 = std::chrono::high_resolution_clock::now();

	#define T  float

	iv_t<T> iv {0.0f, 0.0f, -1.0f};
	timeBound_t<T> tb {0.0f, 2.0f};

	function_t<T> f1([](std::vector<T> f) {
		T t = f.at(0);
		T u1 = f.at(1);
		T u2 = f.at(2);

		return -4 * u1 - 2 * u2 + cosf(t) + 4 * sinf(t);
	});


	function_t<T> f2([](std::vector<T> f) {
		T t = f.at(0);
		T u1 = f.at(1);
		T u2 = f.at(2);

		return 3 * u1 + u2 - 3 * sinf(t);
	});

	ODESystem<T> ode (iv, {f1, f2}, tb, 0.001);

	DataFrame<T> data = solve(ode, ALGORITHM_EULER);

	auto t2 = std::chrono::high_resolution_clock::now();

	std::cout << data << std::endl;

	std::chrono::duration<double, std::milli> ms = t2 - t1;
	std::cout << ms.count() << " ms" << std::endl;

	return 0;
}