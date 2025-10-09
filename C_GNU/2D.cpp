#include <mpfr.h>
#include <vector>
#include <iostream>
using namespace std;


const int PRECISION = 256;
static double EPSILON = 0.01;
static double c = 3.0e8;
static double G = 6.674e-11;

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
			str += "e";
			str += to_string(exp - 1);
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
	
	string ToString(int base = 10, int precision = 10) {
		string str;
		mp_exp_t exp;
		str += this->x.ToString(base, precision) + ", ";
		str += this->y.ToString(base, precision) + ", ";
		str += this->z.ToString(base, precision);
		return "(" + str + ")";
	}
	
	~vec3() {
		try {
			this->x.clear();
			this->y.clear();
			this->z.clear();
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
	std::cout << "M: " << mass.ToString(10, 5) << ", G: " << G << ", r: " << (pos->y).ToString(10, 5) << ", c: " << c << std::endl;
	floating gammaTTR = floating(((G*mass)/((*pos->r)*(*pos->r)*c*c))*(1.0/(1.0 - (2.0*G*mass/((*pos->r)*c*c)))));
	std::cout << "gammaTTR: " << gammaTTR.ToString(10, 10) << std::endl;

	floating gammaTRT = floating(gammaTTR);
	std::cout << "gammaTRT: " << gammaTRT.ToString(10, 10) << std::endl;

	floating gammaRTT = floating((G*mass/((*pos->r)*(*pos->r)))*(1.0 - (2.0*G*mass/((*pos->r)*c*c))));
	std::cout << "gammaRTT: " << gammaRTT.ToString(10, 10) << std::endl;

	floating gammaRRR = floating(-G*mass/((*pos->r)*(*pos->r)*c*c)*(1.0/(1.0 - (2.0*G*mass/((*pos->r)*c*c)))));
	std::cout << "gammaRRR: " << gammaRRR.ToString(10, 10) << std::endl;


	// temp = 1-2GM/(rc^2)
	// mpfr_mul_d(temp, mass, G, MPFR_RNDN);
	// mpfr_mul_d(temp, temp, 2, MPFR_RNDN);
	// mpfr_set_d(temp2, c, MPFR_RNDN);
	// mpfr_sqr(temp2, temp2, MPFR_RNDN);
	// mpfr_mul(temp2, temp2, (*pos->r), MPFR_RNDN);
	// mpfr_div(temp, temp, temp2, MPFR_RNDN);
	// mpfr_d_sub(temp, 1, temp, MPFR_RNDN);

	// gammaRTT = GM/r^2
	// mpfr_mul_d(gammaRTT, mass, G, MPFR_RNDN);
	// mpfr_sqr(temp2, (*pos->r), MPFR_RNDN);
	// mpfr_div(gammaRTT, gammaRTT, temp2, MPFR_RNDN);
	// gammaRTT = GM/r^2(1-2GM/(rc^2))
	// mpfr_mul(gammaRTT, gammaRTT, temp, MPFR_RNDN);


	// temp = 1/(1-2GM/(rc^2))
	// mpfr_d_div(temp, 1, temp, MPFR_RNDN);

	// gammaTTR = GM/(r^2c^2)
	// mpfr_set_d(temp2, c, MPFR_RNDN);
	// mpfr_mul(temp2, temp2, (*pos->r), MPFR_RNDN);
	// mpfr_sqr(temp2, temp2, MPFR_RNDN);
	// mpfr_mul_d(gammaTTR, mass, G, MPFR_RNDN);
	// mpfr_div(gammaTTR, gammaTTR, temp2, MPFR_RNDN);

	// gammaRRR = gammaTRT = gammaTTR = GM/(r^2c^2)(1/(1-2GM/(rc^2)))
	// mpfr_mul(gammaTTR, gammaTTR, temp, MPFR_RNDN);
	// mpfr_set(gammaTRT, gammaTTR, MPFR_RNDN);
	//gammaRRR = -GM/(r^2c^2)(1/(1-2GM/(rc^2)))
	// mpfr_mul_d(gammaRRR, gammaRRR, -1, MPFR_RNDN);
	
	//vec3 acceleration = vec3(0.0, 0.0, 0.0);
	
	// // acceleration.t = (gammaTTR + gammaTRT)vel.t * vel.r
	*acceleration->t = (-gammaTTR * (*vel->t)*(*vel->r)) - (gammaTRT * (*vel->t)*(*vel->r));
	// mpfr_add(acceleration->x, gammaTTR, gammaTRT, MPFR_RNDN);
	// mpfr_mul(acceleration->x, acceleration->x, vel->x, MPFR_RNDN);
	// mpfr_mul(acceleration->x, acceleration->x, vel->y, MPFR_RNDN);
	
	// // acceleration.r = gammaRTT * vel.t^2 - gammaRRR * vel.r^2
	*acceleration->r = (-gammaRTT * (*vel->t)*(*vel->t)) - (gammaRRR * (*vel->r)*(*vel->r));
	// std::cout << "gammaRTT: " << (-gammaRTT * (*vel->t)*(*vel->t)).ToString(10, 10) << " - " << ((gammaRRR * (*vel->r)*(*vel->r))).ToString(10, 10) << " = " << (*acceleration->r).ToString(10, 10) << std::endl;
	// std::cout << "Acceleration: " << acceleration->ToString(10, 10) << std::endl;
	// mpfr_sqr(temp, vel->x , MPFR_RNDN);
	// mpfr_mul(temp, temp, gammaRTT, MPFR_RNDN);
	// mpfr_sqr(acceleration->y, vel->y, MPFR_RNDN);
	// mpfr_mul(acceleration->y, acceleration->y, gammaRRR, MPFR_RNDN);
	// mpfr_sub(acceleration->y, temp, acceleration->y, MPFR_RNDN);

	gammaTTR.clear();
	gammaTRT.clear();
	gammaRRR.clear();
	gammaRTT.clear();
	
	// std::cout << vel->y.ToString(10, 10) << " + "<< acceleration->y.ToString(10, 10) << " = ";
	// std::cout << "acceleration: " << acceleration->y.ToString(10, 10) << std::endl;
	(*vel) += (*acceleration);
	std::cout << vel->y.ToString(10, 10) << std::endl;

	(*pos) += (*vel);
	std::cout << "pos: " << pos->ToString() << std::endl;
}

int main() {
	vec3 pos = vec3(0, 10, 0);
	std::cout << "pos orig: " << pos.ToString() << std::endl;

	vec3 vel = vec3(EPSILON, 0, 0);

	floating mass = floating(100.0, PRECISION);

	vec3 acceleration = vec3(0.0, 0.0, 0.0);

	for (int i = 0; i < 20; i++) {
		step(&pos, &vel, mass, &acceleration);
	}

	// Exit
	return 0;
}