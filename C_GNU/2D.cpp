#include <mpfr.h>
#include <vector>
#include <iostream>
using namespace std;


const int PRECISION = 256;
const double EPSILON = 0.01;
const double c = 3.0e8;
const double G = 6.674e-11;

class floating {
	public:
	mpfr_t num;
	floating() {
		mpfr_init2(this->num, PRECISION);
		mpfr_set_d(this->num, 0.0, MPFR_RNDN);
	}
	floating(double num, int precision = PRECISION) {
		mpfr_init2(this->num, precision);
		mpfr_set_d(this->num, num, MPFR_RNDN);
	}
	floating(floating& num, int precision = 0) {
		precision = precision == 0 ? mpfr_get_prec(num.num) : precision;
		mpfr_init2(this->num, precision);
		mpfr_set(this->num, num.num, MPFR_RNDN);
	}
	floating(mpfr_t num, int precision = 0) {
		precision = precision == 0 ? mpfr_get_prec(num) : precision;
		mpfr_init2(this->num, precision);
		mpfr_set(this->num, num, MPFR_RNDN);
	}
	void clear() {
		try {
			mpfr_clear(this->num);
		} catch (...) {
			std::cout << "floating already deallocated";
		}
	}



	operator string() {
		return ToString();
	}

	string ToString(int base = 10, int precision = 0) {
		// base: 2–62 (10 for decimal)
		// precision: number of digits; 0 = use full precision

		mp_exp_t exp;
		char *mantissa = mpfr_get_str(nullptr, &exp, base, precision, this->num, MPFR_RNDN);

		if (!mantissa)
			return "<conversion error>";

		string str;
		//str += mantissa[0];
		int i = 0;
		if (mpfr_sgn(this->num) < 0) {
			str += "-";
			i++;
		}

		// Build number: mantissa + exponent
		if (mantissa[i] != '\0') {
			str += (mantissa[i]);
			str += '.';
			str += (mantissa + i + 1);
		}
		// Want all negative exponents to show, but not 0 or 1
		if (exp != 1) {
			str += "*10^{";
			str += to_string(exp - 1) + "}";
		}

		mpfr_free_str(mantissa);
		return str;
	}
	
	void operator =(mpfr_t other) {
		mpfr_set(this->num, other, MPFR_RNDN);
	}
	void operator =(double other) {
		mpfr_set_d(this->num, other, MPFR_RNDN);
	}

	floating operator +(floating other) {
		floating ret = floating(0.0);
		mpfr_add(ret.num, this->num, other.num, MPFR_RNDN);
		return ret;
	}
	floating operator +(double other) {
		floating ret = floating(0.0);
		mpfr_add_d(ret.num, this->num, other, MPFR_RNDN);
		return ret;
	}
	void operator +=(floating other) {
		mpfr_add(this->num, this->num, other.num, MPFR_RNDN);
	}

	
	floating operator -(floating other) {
		floating ret = floating(0.0);
		mpfr_sub(ret.num, this->num, other.num, MPFR_RNDN);
		return ret;
	}
	floating operator -(double other) {
		floating ret = floating(0.0);
		mpfr_sub_d(ret.num, this->num, other, MPFR_RNDN);
		return ret;
	}


	floating operator *(floating other) {
		floating ret = floating(0.0);
		mpfr_mul(ret.num, this->num, other.num, MPFR_RNDN);
		return ret;
	}
	floating operator *(double other) {
		floating ret = floating(0.0);
		mpfr_mul_d(ret.num, this->num, other, MPFR_RNDN);
		return ret;
	}
	void operator *=(double other) {
		mpfr_mul_d(this->num, this->num, other, MPFR_RNDN);
	}


	floating operator /(floating other) {
		floating ret = floating(0.0);
		mpfr_div(ret.num, this->num, other.num, MPFR_RNDN);
		return ret;
	}
	floating operator /(double other) {
		floating ret = floating(0.0);
		mpfr_div_d(ret.num, this->num, other, MPFR_RNDN);
		return ret;
	}
};

floating operator -(double left, floating right) {
	floating ret = floating(0.0, mpfr_get_prec(right.num));
	mpfr_d_sub(ret.num, left, right.num, MPFR_RNDN);
	return ret;
}
floating operator -(floating operand) {
	return operand * -1.0;
}

floating operator *(double left, floating right) {
	floating ret = floating(0.0, mpfr_get_prec(right.num));
	mpfr_mul_d(ret.num, right.num, left, MPFR_RNDN);
	return ret;
}

floating operator /(double left, floating right) {
	floating ret = floating(0.0, mpfr_get_prec(right.num));
	mpfr_d_div(ret.num, left, right.num, MPFR_RNDN);
	return ret;
}

class vec3 {
public:
	string name;
	// (x,y,z), (t,r,phi)
	floating x;
	floating y;
	floating z;

	floating* t;
	floating* r;
	floating* phi;

	vec3(double xt, double yr, double z, int precision = PRECISION) {
		this->x = floating(xt, precision);
		this->y = floating(yr, precision);
		this->z = floating(z, precision);

		this->t = &this->x;
		this->r = &this->y;
		this->phi = &this->z;
	}
	vec3(floating xt, floating yr, floating z, int precision = PRECISION) {
		this->x = xt;
		this->y = yr;
		this->z = z;

		this->t = &this->x;
		this->r = &this->y;
		this->phi = &this->z;
	}
	vec3(const vec3& vec) {
		this->x = floating(vec.x);
		this->y = floating(vec.y);
		this->z = floating(vec.z);

		this->t = &this->x;
		this->r = &this->y;
		this->phi = &this->z;
	}
	
	string ToString(int base = 10, int precision = 10) {
		string str;
		mp_exp_t exp;
		str += this->x.ToString(base, precision) + ", ";
		str += this->y.ToString(base, precision);
		return "(" + str + ")";
	}
	
	~vec3() {
		try {
			x.clear();
			y.clear();
			z.clear();
			t = nullptr;
			r = nullptr;
			phi = nullptr;
		} catch (std::exception& ex) {
			std::cout << ": vec3 already deallocated" << std::endl;
		} catch (std::string& ex) {
			std::cout << ": vec3 already deallocated" << std::endl;	
		} catch (...) {
			std::cout << ": vec3 already deallocated" << std::endl;
		}
	}

	vec3 operator +(vec3 other) {
		return vec3(this->x + other.x, this->y + other.y, this->z + other.z);
	}
	void operator +=(vec3 other) {
		this->x += other.x;
		this->y += other.y;
		this->z += other.z;
	}
	vec3 operator -(vec3 other) {
		return vec3(this->x - other.x, this->y - other.y, this->z - other.z);
	}
	vec3 operator *(floating other) {
		return vec3(this->x * other, this->y * other, this->z * other);
	}
	vec3 operator /(floating other) {
		return vec3(this->x / other, this->y + other, this->z + other);
	}
};



static void step(vec3* pos, vec3* vel, floating mass, vec3* acceleration) {


	// Calculate Christoffel Symbols

	// GM/(r^2c^2)(1/(1-2GM/(rc^2)))
	// std::cout << "M: " << mass.ToString(10, 5) << ", G: " << G << ", r: " << (pos->y).ToString(10, 5) << ", c: " << c << std::endl;
	floating gammaTTR = floating(((G*mass)/((*pos->r)*(*pos->r)*c*c))*(1.0/(1.0 - (2.0*G*mass/((*pos->r)*c*c)))));
	// std::cout << "gammaTTR: " << gammaTTR.ToString(10, 10) << std::endl;

	floating gammaTRT = floating(gammaTTR);
	// std::cout << "gammaTRT: " << gammaTRT.ToString(10, 10) << std::endl;

	floating gammaRTT = floating((G*mass/((*pos->r)*(*pos->r)))*(1.0 - (2.0*G*mass/((*pos->r)*c*c))));
	// std::cout << "gammaRTT: " << gammaRTT.ToString(10, 10) << std::endl;

	floating gammaRRR = floating(-G*mass/((*pos->r)*(*pos->r)*c*c)*(1.0/(1.0 - (2.0*G*mass/((*pos->r)*c*c)))));
	// std::cout << "gammaRRR: " << gammaRRR.ToString(10, 10) << std::endl;

	
	// // acceleration.t = (gammaTTR + gammaTRT)vel.t * vel.r
	*acceleration->t = (-gammaTTR * (*vel->t)*(*vel->r)) - (gammaTRT * (*vel->t)*(*vel->r));
	
	// // acceleration.r = gammaRTT * vel.t^2 - gammaRRR * vel.r^2
	*acceleration->r = (-gammaRTT * (*vel->t)*(*vel->t)) - (gammaRRR * (*vel->r)*(*vel->r));
	// std::cout << "gammaRTT: " << (-gammaRTT * (*vel->t)*(*vel->t)).ToString(10, 10) << " - " << ((gammaRRR * (*vel->r)*(*vel->r))).ToString(10, 10) << " = " << (*acceleration->r).ToString(10, 10) << std::endl;
	// std::cout << "Acceleration: " << acceleration->ToString(10, 10) << std::endl;


	gammaTTR.clear();
	gammaTRT.clear();
	gammaRRR.clear();
	gammaRTT.clear();
	
	// std::cout << vel->y.ToString(10, 10) << " + "<< acceleration->y.ToString(10, 10) << " = ";
	// std::cout << "acceleration: " << acceleration->y.ToString(10, 10) << std::endl;
	(*vel) += (*acceleration);
	// std::cout << "Vel: " << vel->ToString(10, 10) << std::endl;
	// std::cout << vel->y.ToString(10, 10) << std::endl;

	(*pos) += (*vel);
	std::cout << pos->ToString(10, 10) << ", ";
}

int main() {
	vec3 pos = vec3(0, 2, 0);
	std::cout << pos.ToString() << ", ";

	vec3 vel = vec3(EPSILON, 0.001, 0);

	floating mass = floating(1000000000.0);

	vec3 acceleration = vec3(0.0, 0.0, 0.0);

	for (int i = 0; i < 10000; i++) {
		step(&pos, &vel, mass, &acceleration);
		if (mpfr_get_d(pos.y.num, MPFR_RNDN) <= 0.0) {
			break;
		}
	}

	// Exit
	return 0;
}