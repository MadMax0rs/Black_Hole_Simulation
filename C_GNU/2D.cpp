#include <mpfr.h>
#include <vector>
#include <iostream>
using namespace std;


const int PRECISION = 256;
static double EPSILON = 0.0001;
static double c = 3.0e8;
static double G = 6.674e-11;

class vec3 {
public:
	string name;
	// (x,y,z), (t,r,phi)
	mpfr_t x;
	mpfr_t y;
	mpfr_t z;
	bool initialized;

	vec3(double x, double y, double z, int precision = PRECISION) {
		mpfr_init2(this->x, precision);
		mpfr_init2(this->y, precision);
		mpfr_init2(this->z, precision);

		mpfr_set_d(this->x, x, MPFR_RNDN);
		mpfr_set_d(this->y, y, MPFR_RNDN);
		mpfr_set_d(this->z, z, MPFR_RNDN);
		this->initialized = true;
	}
	vec3(mpfr_t x, mpfr_t y, mpfr_t z, int precision = PRECISION) {
		mpfr_init2(this->x, precision);
		mpfr_init2(this->y, precision);
		mpfr_init2(this->z, precision);
		
		mpfr_set(this->x, x, MPFR_RNDN);
		mpfr_set(this->y, y, MPFR_RNDN);
		mpfr_set(this->z, z, MPFR_RNDN);
		this->initialized = true;
	}
	
	void print() {
		char* str;
		mp_exp_t exp;
		str = mpfr_get_str(NULL, &exp, 10, 20, x, MPFR_RNDN);
		std::cout << "(" << str << " * 10^" << exp << ")";
		mpfr_free_str(str);
	}
	
	~vec3() {
		try {
			if (!this->initialized) {
				return;
			}
			mpfr_clear(this->x);
			mpfr_clear(this->y);
			mpfr_clear(this->z);
			this->initialized = false;
		} catch (std::exception& ex) {
			std::cout << ": vec3 already initialized" << std::endl;
		} catch (std::string& ex) {
			std::cout << ": vec3 already initialized" << std::endl;	
		} catch (...) {
			std::cout << ": vec3 already initialized" << std::endl;
		}
	}

	vec3 operator +(vec3 other) {
		vec3 ret = vec3(0.0, 0.0, 0.0);
		mpfr_add(ret.x, x, other.x, MPFR_RNDN);
		mpfr_add(ret.y, y, other.y, MPFR_RNDN);
		mpfr_add(ret.z, z, other.z, MPFR_RNDN);
		return ret;
	}
	void operator +=(vec3 other) {
		mpfr_add(x, x, other.x, MPFR_RNDN);
		mpfr_add(y, y, other.y, MPFR_RNDN);
		mpfr_add(z, z, other.z, MPFR_RNDN);
	}
	vec3 operator -(vec3 other) {
		vec3 ret = vec3(0.0, 0.0, 0.0);
		mpfr_sub(ret.x, x, other.x, MPFR_RNDN);
		mpfr_sub(ret.y, y, other.y, MPFR_RNDN);
		mpfr_sub(ret.z, z, other.z, MPFR_RNDN);
		return ret;
	}
	vec3 operator *(mpfr_t other) {
		vec3 ret = vec3(0.0, 0.0, 0.0);
		mpfr_mul(ret.x, x, other, MPFR_RNDN);
		mpfr_mul(ret.y, y, other, MPFR_RNDN);
		mpfr_mul(ret.z, z, other, MPFR_RNDN);
		return ret;
	}
	vec3 operator /(mpfr_t other) {
		vec3 ret = vec3(0.0, 0.0, 0.0);
		mpfr_div(ret.x, x, other, MPFR_RNDN);
		mpfr_div(ret.y, y, other, MPFR_RNDN);
		mpfr_div(ret.z, z, other, MPFR_RNDN);
		return ret;
	}
};



static void step(vec3* pos, vec3* vel, mpfr_t mass, vec3* acceleration) {


	// Calculate Christoffel Symbols
	mpfr_t temp;
	mpfr_init2(temp, PRECISION);
	mpfr_t temp2;
	mpfr_init2(temp2, PRECISION);
	mpfr_t gammaTTR;
	mpfr_init2(gammaTTR, PRECISION);
	mpfr_t gammaTRT;
	mpfr_init2(gammaTRT, 64);
	mpfr_t gammaRTT;
	mpfr_init2(gammaRTT, PRECISION);
	mpfr_t gammaRRR;
	mpfr_init2(gammaRRR, PRECISION);


	// temp = 1-2GM/(rc^2)
	mpfr_mul_d(temp, mass, G, MPFR_RNDN);
	mpfr_mul_d(temp, temp, 2, MPFR_RNDN);
	mpfr_set_d(temp2, c, MPFR_RNDN);
	mpfr_sqr(temp2, temp2, MPFR_RNDN);
	mpfr_mul(temp2, temp2, pos->y, MPFR_RNDN);
	mpfr_div(temp, temp, temp2, MPFR_RNDN);
	mpfr_d_sub(temp, 1, temp, MPFR_RNDN);

	// gammaRTT = GM/r^2
	mpfr_mul_d(gammaRTT, mass, G, MPFR_RNDN);
	mpfr_sqr(temp2, pos->y, MPFR_RNDN);
	mpfr_div(gammaRTT, gammaRTT, temp2, MPFR_RNDN);
	// gammaRTT = GM/r^2(1-2GM/(rc^2))
	mpfr_mul(gammaRTT, gammaRTT, temp, MPFR_RNDN);


	// temp = 1/(1-2GM/(rc^2))
	mpfr_d_div(temp, 1, temp, MPFR_RNDN);

	// gammaTTR = GM/(r^2c^2)
	mpfr_set_d(temp2, c, MPFR_RNDN);
	mpfr_mul(temp2, temp2, pos->y, MPFR_RNDN);
	mpfr_sqr(temp2, temp2, MPFR_RNDN);
	mpfr_mul_d(gammaTTR, mass, G, MPFR_RNDN);
	mpfr_div(gammaTTR, gammaTTR, temp2, MPFR_RNDN);

	// // gammaRRR = gammaTRT = gammaTTR = GM/(r^2c^2)(1/(1-2GM/(rc^2)))
	mpfr_mul(gammaTTR, gammaTTR, temp, MPFR_RNDN);
	std::cout << "1";
	mpfr_set(gammaTRT, gammaTTR, MPFR_RNDN);
	// //gammaRRR = -GM/(r^2c^2)(1/(1-2GM/(rc^2)))
	mpfr_mul_d(gammaRRR, gammaRRR, -1, MPFR_RNDN);
	
	//vec3 acceleration = vec3(0.0, 0.0, 0.0);
	
	// acceleration.x = (gammaTTR + gammaTRT)vel.t * vel.r
	mpfr_add(acceleration->x, gammaTTR, gammaTRT, MPFR_RNDN);
	mpfr_mul(acceleration->x, acceleration->x, vel->x, MPFR_RNDN);
	mpfr_mul(acceleration->x, acceleration->x, vel->y, MPFR_RNDN);
	
	// acceleration.y = gammaRTT * vel.t^2 - gammaRRR * vel.r^2
	mpfr_sqr(temp, vel->x , MPFR_RNDN);
	mpfr_mul(temp, temp, gammaRTT, MPFR_RNDN);
	mpfr_sqr(acceleration->y, vel->y, MPFR_RNDN);
	mpfr_mul(acceleration->y, acceleration->y, gammaRRR, MPFR_RNDN);
	mpfr_sub(acceleration->y, temp, acceleration->y, MPFR_RNDN);
	std::cout << "2";

	mpfr_clear(gammaTTR);
	mpfr_clear(gammaTRT);
	mpfr_clear(gammaRRR);
	mpfr_clear(gammaRTT);
	mpfr_clear(temp);
	mpfr_clear(temp2);
}

int main() {
	vec3 pos = vec3(0, -20, 0);
	pos.name = "pos";
	vec3 vel = vec3(0, c, 0);
	vel.name = "vel";
	mpfr_t mass;
	mpfr_init2(mass, PRECISION);
	mpfr_set_d(mass, 1.0, MPFR_RNDN);
	vec3 acceleration = vec3(0.0, 0.0, 0.0);
	acceleration.name = "acceleration";
	step(&pos, &vel, mass, &acceleration);

	vel += acceleration;
	pos += vel;

	
	pos.print();

	// Exit
	return 0;
}