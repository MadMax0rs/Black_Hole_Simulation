#include <mpfr.h>
#include <iostream>

using namespace std;

const int PRECISION = 256;

string mpfrToString(const mpfr_t x, int base = 10, int precision = 0) {
    // base: 2–62 (10 for decimal)
    // precision: number of digits; 0 = use full precision

    mp_exp_t exp;
    char *mantissa = mpfr_get_str(nullptr, &exp, base, precision, x, MPFR_RNDN);

    if (!mantissa)
        return "<conversion error>";

    string str;
	str += mantissa[0];
    //if (mpfr_sgn(x) < 0) str += "-";

    // Build number: mantissa + exponent
	if (mantissa[1] != '\0') {
		str += (mantissa[1]);
		str += '.';
		str += (mantissa + 2);
	}
	str += "e";
	str += to_string(exp - 1);

    mpfr_free_str(mantissa);
    return str;
}


int main() {
	mpfr_t number;
	mpfr_init2(number, PRECISION);
	mpfr_set_d(number, -0.2050505051, MPFR_RNDN);

	std::cout << mpfrToString(number, 10, 10);
	return 0;
}