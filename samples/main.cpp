// #include <glad/glad.h>
// #include <GLFW/glfw3.h>

#include <iostream>
#include <cmath>
#include "diffeq.h"


using namespace DES;


int main() {

	std::function<float(std::vector<float>)> func1 = [](std::vector<float> f) {
		float t = f.at(0);
		float u1 = f.at(1);
		float u2 = f.at(2);

		return -4 * u1 - 2 * u2 + cos(t) + 4 * sin(t);
	};

	std::function<float(std::vector<float>)> func2 = [](std::vector<float> f) {
		float t = f.at(0);
		float u1 = f.at(1);
		float u2 = f.at(2);

		return 3 * u1 + u2 - 3 * sin(t);
	};

	iv_t iv {0.0f, 0.0f, -1.0f};
	timeBound_t tb {0.0f, 0.2f};

	function_t<float, float> f1(func1);
	function_t<float, float> f2(func2);

	ODESystem<float, float> ode (iv, {f1, f2}, tb, 0.1);

	auto data = solve<float, float>(ode, 3);


	return 0;
}