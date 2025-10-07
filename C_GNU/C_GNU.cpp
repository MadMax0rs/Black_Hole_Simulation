#include <mpfr.h>
#include <vector>
using namespace std;


const int PRECISION = 256;
static mpfr_t EPSILON;
static mpfr_t c;
static mpfr_t G;

class vec3 {
	public:
		mpfr_t x;
		mpfr_t y;
		mpfr_t z;
		
		vec3(double x, double y, double z, int precision = PRECISION) {
			mpfr_init2(this->x, precision);
			mpfr_init2(this->y, precision);
			mpfr_init2(this->z, precision);

			mpfr_set_d(this->x, x, MPFR_RNDN);
			mpfr_set_d(this->y, y, MPFR_RNDN);
			mpfr_set_d(this->z, z, MPFR_RNDN);
		}
		vec3(mpfr_t x, mpfr_t y, mpfr_t z, int precision = PRECISION) {
			mpfr_init2(this->x, precision);
			mpfr_init2(this->y, precision);
			mpfr_init2(this->z, precision);

			mpfr_set(this->x, x, MPFR_RNDN);
			mpfr_set(this->y, y, MPFR_RNDN);
			mpfr_set(this->z, z, MPFR_RNDN);
		}

		~vec3() {
			mpfr_clear(x);
			mpfr_clear(y);
			mpfr_clear(z);
		}

	vec3 operator +(vec3 other) {
		vec3 ret = vec3(0.0,0.0,0.0);
		mpfr_add(ret.x, x, other.x, MPFR_RNDN);
		mpfr_add(ret.y, y, other.y, MPFR_RNDN);
		mpfr_add(ret.z, z, other.z, MPFR_RNDN);
		return ret;
	}
	vec3 operator -(vec3 other) {
		vec3 ret = vec3(0.0,0.0,0.0);
		mpfr_sub(ret.x, x, other.x, MPFR_RNDN);
		mpfr_sub(ret.y, y, other.y, MPFR_RNDN);
		mpfr_sub(ret.z, z, other.z, MPFR_RNDN);
		return ret;
	}
	vec3 operator *(mpfr_t other) {
		vec3 ret = vec3(0.0,0.0,0.0);
		mpfr_mul(ret.x, x, other, MPFR_RNDN);
		mpfr_mul(ret.y, y, other, MPFR_RNDN);
		mpfr_mul(ret.z, z, other, MPFR_RNDN);
		return ret;
	}
	vec3 operator /(mpfr_t other) {
		vec3 ret = vec3(0.0,0.0,0.0);
		mpfr_div(ret.x, x, other, MPFR_RNDN);
		mpfr_div(ret.y, y, other, MPFR_RNDN);
		mpfr_div(ret.z, z, other, MPFR_RNDN);
		return ret;
	}
};

class vec4 {
	public:
		mpfr_t t;
		mpfr_t r;
		mpfr_t theta;
		mpfr_t phi;
		
		vec4(double t, double r, double theta, double phi, int precision = 256) {
			mpfr_init2(this->t, precision);
			mpfr_init2(this->r, precision);
			mpfr_init2(this->theta, precision);
			mpfr_init2(this->phi, precision);

			mpfr_set_d(this->t, t, MPFR_RNDN);
			mpfr_set_d(this->r, r, MPFR_RNDN);
			mpfr_set_d(this->theta, theta, MPFR_RNDN);
			mpfr_set_d(this->phi, phi, MPFR_RNDN);
		}
		vec4(mpfr_t t, mpfr_t r, mpfr_t theta, mpfr_t phi, int precision = 256) {
			mpfr_init2(this->t, precision);
			mpfr_init2(this->r, precision);
			mpfr_init2(this->theta, precision);
			mpfr_init2(this->phi, precision);

			mpfr_set(this->t, t, MPFR_RNDN);
			mpfr_set(this->r, r, MPFR_RNDN);
			mpfr_set(this->theta, theta, MPFR_RNDN);
			mpfr_set(this->phi, phi, MPFR_RNDN);
		}

		// Copy constructor
		vec4(const vec4& other) {
			mpfr_init2(t, mpfr_get_prec(other.t)); mpfr_set(t, other.t, MPFR_RNDN);
			mpfr_init2(r, mpfr_get_prec(other.r)); mpfr_set(r, other.r, MPFR_RNDN);
			mpfr_init2(theta, mpfr_get_prec(other.theta)); mpfr_set(theta, other.theta, MPFR_RNDN);
			mpfr_init2(phi, mpfr_get_prec(other.phi)); mpfr_set(phi, other.phi, MPFR_RNDN);
		}

		~vec4() {
			mpfr_clear(t);
			mpfr_clear(r);
			mpfr_clear(theta);
			mpfr_clear(phi);
		}

	vec4 operator +(vec4 other) {
		vec4 ret = vec4(0.0,0.0,0.0,0.0);
		mpfr_add(ret.t, t, other.t, MPFR_RNDN);
		mpfr_add(ret.r, r, other.r, MPFR_RNDN);
		mpfr_add(ret.theta, theta, other.theta, MPFR_RNDN);
		mpfr_add(ret.phi, phi, other.phi, MPFR_RNDN);
		return ret;
	}
	vec4 operator +(mpfr_t other) {
		vec4 ret = vec4(0.0,0.0,0.0,0.0);
		mpfr_add(ret.t, t, other, MPFR_RNDN);
		mpfr_add(ret.r, r, other, MPFR_RNDN);
		mpfr_add(ret.theta, theta, other, MPFR_RNDN);
		mpfr_add(ret.phi, phi, other, MPFR_RNDN);
		return ret;
	}
	vec4 operator -(vec4 other) {
		vec4 ret = vec4(0.0,0.0,0.0,0.0);
		mpfr_sub(ret.t, t, other.t, MPFR_RNDN);
		mpfr_sub(ret.r, r, other.r, MPFR_RNDN);
		mpfr_sub(ret.theta, theta, other.theta, MPFR_RNDN);
		mpfr_sub(ret.phi, phi, other.phi, MPFR_RNDN);
		return ret;
	}
	vec4 operator *(mpfr_t other) {
		vec4 ret = vec4(0.0,0.0,0.0,0.0);
		mpfr_mul(ret.t, t, other, MPFR_RNDN);
		mpfr_mul(ret.r, r, other, MPFR_RNDN);
		mpfr_mul(ret.theta, theta, other, MPFR_RNDN);
		mpfr_mul(ret.phi, phi, other, MPFR_RNDN);
		return ret;
	}
	vec4 operator /(mpfr_t other) {
		vec4 ret = vec4(0.0,0.0,0.0,0.0);
		mpfr_div(ret.t, t, other, MPFR_RNDN);
		mpfr_div(ret.r, r, other, MPFR_RNDN);
		mpfr_div(ret.theta, theta, other, MPFR_RNDN);
		mpfr_div(ret.phi, phi, other, MPFR_RNDN);
		return ret;
	}
	vec3 yzw() {
		return vec3(r, theta, phi);
	}
};

class State {
	public:
		vec4 rk4;
		vec4 p;

		State(vec4 rk4, vec4 p) :
			rk4(rk4),
			p(p)
		{}
	
	State operator +(mpfr_t other) {
		return State(rk4+other, p+other);
	}
	State operator *(mpfr_t other) {
		return State(rk4*other, p*other);
	}

	vector<vec4> split() const {
		return {rk4, p};
	}
};

struct Ray {
public:
	vec3 pos;
	vec3 col;
	int numSteps;

	Ray(vec3 pos, vec3 col, int numSteps) :
		pos(pos),
		col(col),
		numSteps(numSteps)
	{}
};

int main() {
	mpfr_init2(EPSILON, PRECISION);
	mpfr_init2(c, PRECISION);
	mpfr_init2(G, PRECISION);
	mpfr_set_d(EPSILON, 0.0001, MPFR_RNDN);
	mpfr_set_d(c, 3.0e8, MPFR_RNDN);
	mpfr_set_d(G, 6.674e-11, MPFR_RNDN);

	vec3 vec = vec3(1, 2, 3);
	Ray ray = Ray(vec3(0, 0, -20), vec3(0.0, 0.0, 0.0), 0);
	


	mpfr_t gammaTTR;
	mpfr_t gammaTRT;
	mpfr_t gammaRTT;
	mpfr_t gammaRRR;
	mpfr_t gammaRThetaTheta;
	mpfr_t gammaRPhiPhi;
	mpfr_t gammaThetaRTheta;
	mpfr_t gammaThetaThetaR;
	mpfr_t gammaThetaPhiPhi;
	mpfr_t gammaPhiRPhi;
	mpfr_t gammaPhiPhiR;
	mpfr_t gammaPhiThetaPhi;
	mpfr_t gammaPhiPhiTheta;


	// Exit
    return 0;
}