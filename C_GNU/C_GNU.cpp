#include <mpfr.h>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;

class vec4;

const int PRECISION = 64;
const double EPSILON = 0.001;
const double c = 1; //3.0e8;
const double G = 1; //1.7518e-45; //6.674e-11;


class floating {
public:
	mpfr_t num;
	string name;
	bool deallocated;
	floating() {
		mpfr_init2(this->num, PRECISION);
		mpfr_set_d(this->num, 0.0, MPFR_RNDN);
		deallocated = false;
	}
	floating(double num, int precision = PRECISION) {
		mpfr_init2(this->num, precision);
		mpfr_set_d(this->num, num, MPFR_RNDN);
		deallocated = false;
	}
	floating(const floating& num, int precision = 0) {
		precision = precision == 0 ? mpfr_get_prec(num.num) : precision;
		mpfr_init2(this->num, precision);
		mpfr_set(this->num, num.num, MPFR_RNDN);
		deallocated = false;
	}
	floating(mpfr_t num, int precision = 0) {
		precision = precision == 0 ? mpfr_get_prec(num) : precision;
		mpfr_init2(this->num, precision);
		mpfr_set(this->num, num, MPFR_RNDN);
		deallocated = false;
	}
	void clear() {
		try {
			if (deallocated) {
				std::cerr << name << " was already deallocated" << std::endl;
				return;
			}
			mpfr_clear(this->num);
			deallocated = true;
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
	void operator =(const floating& other) {
		mpfr_init2(num, mpfr_get_prec(other.num));
		mpfr_set(num, other.num, MPFR_RNDN);
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

	bool operator <=(floating other) {
		return mpfr_lessequal_p(this->num, other.num) != 0;
	}
};

	static floating operator /(const floating& left, floating other) {
		floating ret = floating(0.0);
		mpfr_div(ret.num, left.num, other.num, MPFR_RNDN);
		return ret;
	}

	static floating operator /(const floating& left, double other) {
		floating ret = floating(0.0);
		mpfr_div_d(ret.num, left.num, other, MPFR_RNDN);
		return ret;
	}

static floating operator -(double left, floating right) {
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

// Floating trig
const floating sin(floating a) {
	// std::cout << "sin(" << a.ToString() << std::endl;
	floating ret = floating();
	mpfr_sin(ret.num, a.num, MPFR_RNDN);
	return ret;
}
const floating asin(floating a) {
	floating ret = floating();
	mpfr_asin(ret.num, a.num, MPFR_RNDN);
	return ret;
}
const floating cos(const floating& a) {
	floating ret = floating();
	mpfr_cos(ret.num, a.num, MPFR_RNDN);
	return ret;
}
const floating acos(floating a) {
	floating ret = floating();
	mpfr_acos(ret.num, a.num, MPFR_RNDN);
	return ret;
}
const floating atan2(floating opposite, floating adjacent) {
	floating ret = floating();
	mpfr_atan2(ret.num, opposite.num, adjacent.num, MPFR_RNDN);
	return ret;
}
const floating cot(floating a) {
	floating ret = floating();
	mpfr_cot(ret.num, a.num, MPFR_RNDN);
	return ret;
}
const floating sqrt(floating a) {
	floating ret = floating();
	mpfr_sqrt(ret.num, a.num, MPFR_RNDN);
	return ret;
}



class vec3 {
public:
	floating x;
	floating y;
	floating z;

	vec3(double x, double y, double z, int precision = PRECISION) {
		this->x = floating(x, precision);
		this->y = floating(y, precision);
		this->z = floating(z, precision);
	}
	vec3(floating x, floating y, floating z, int precision = PRECISION) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	vec3(const vec4& vec, int precision);
	
	~vec3() {
		clear();
	}

	void clear() {
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
	
	string ToString(int base = 10, int precision = 10, bool xy = true) {
		string str;
		mp_exp_t exp;
		str += this->x.ToString(base, precision) + ", ";
		if (!xy) str += this->y.ToString(base, precision) + ", ";
		str += this->z.ToString(base, precision);
		return "(" + str + ")";
	}

	vec3 toCartesian() {
		vec3 ret = vec3(0,0,0);
		ret.x = x * cos(z) * sin(y);
		ret.y = x * cos(y);
		ret.z = x * sin(z) * sin(y);
		// x = r, y = theta, z = phi
		return ret;
	}


	vec3 toSphere() {
		// x = r, y = theta, z = phi
		floating len = length();
		return vec3(len, acos(y/len), atan2(z, x));
	}

	floating length() {
		return sqrt((x*x) + (y*y) + (z*z));
	}

	vec3 operator +(vec3 other) {
		return vec3(this->x + other.x, this->y + other.y, this->z + other.z);
	}
	void operator +=(vec3 other) {
		x += other.x;
		y += other.y;
		z += other.z;
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

class vec4 {
public:
	// r, theta, phi
	// distance, azemuth (angle on the xy plane), polar angle(elevation)
	floating t;
	floating r;
	floating theta;
	floating phi;
	
	vec4() {
		vec4(0.0, 0.0, 0.0, 0.0);
	}
	vec4(double t, double r, double theta, double phi, int precision = PRECISION) {
	this->t = floating(t, precision);
	this->r = floating(r, precision);
	this->theta = floating(theta, precision);
	this->phi = floating(phi, precision);
	}
	vec4(floating t, floating r, floating theta, floating phi, int precision = PRECISION) {
		this->t = t;
		this->r = r;
		this->theta = theta;
		this->phi = phi;
	}
	vec4(floating t, const vec3& rThetaPhi, int precision = PRECISION) {
		this->t = t;
		this->r = rThetaPhi.x;
		this->theta = rThetaPhi.y;
		this->phi = rThetaPhi.z;
	}

	// Copy constructor
	vec4(vec4& vec) {
		this->t = floating(vec.t);
		this->r = floating(vec.r);
		this->theta = floating(vec.theta);
		this->phi = floating(vec.phi);
	}

	~vec4() {
		clear();
	}

	void clear() {
		t.clear();
		r.clear();
		theta.clear();
		phi.clear();
	}

	void setName(string name) {
		t.name = name + ".t";
		r.name = name + ".r";
		theta.name = name + ".theta";
		phi.name = name + ".phi";
	}

	vec4 operator +(vec4 other) {
		return vec4(t + other.t, r + other.r, theta + other.theta, phi + other.phi);
	}
	void operator +=(vec4 other) {
		t += other.t;
		r += other.r;
		theta += other.theta;
		phi += other.phi;
	}
	vec4 operator +(floating other) {
		return vec4(t + other, r + other, theta + other, phi + other);
	}
	vec4 operator -(vec4 other) {
		return vec4(t - other.t, r - other.r, theta - other.theta, phi - other.phi);
	}
	vec4 operator *(floating other) {
		return vec4(t * other, r * other, theta * other, phi * other);
	}
	vec4 operator /(floating other) {
		return vec4(t / other, r / other, theta / other, phi / other);
	}
	operator vec3() {
		return vec3(r, theta, phi);
	}
};

// C++ requires types to be defined BEFORE they're used, so have to
// define vec4 before making this constructor
vec3::vec3(const vec4& vec, int precision = PRECISION) {
	this->x = vec.r;
	this->y = vec.theta;
	this->z = vec.phi;
}

static void step(vec4* pos, vec4* vel, vec4* acceleration, floating mass, ofstream* outputFile, bool debug, bool debug2) {
	floating common = pos->r - (2.0 * G * mass);
	// std::cout << "pos.t: " << pos->t.ToString(10, 5) << std::endl;
	std::cout << "pos.r: " << pos->r.ToString(10, 5) << std::endl;
	// std::cout << "pos.theta: " << pos->theta.ToString(10, 5) << std::endl;
	// std::cout << "pos.phi: " << pos->phi.ToString(10, 5) << std::endl;
	// std::cout << "common: " << common.ToString(10, 5) << std::endl;

	floating gammaTTR			= G*mass / (pos->r * common);
	floating gammaTRT			= floating(gammaTTR);

	floating gammaRTT			= G*mass*common / (pos->r * pos->r * pos->r);
	floating gammaRRR			= -gammaTTR;
	floating gammaRThetaTheta	= -common;
	floating gammaRPhiPhi		= gammaRThetaTheta*sin(pos->theta)*sin(pos->theta);
	// std::cout << "testttttt: " << std::endl;
	// std::cout << "\tvel->theta: " << vel->theta.ToString(10, 10) << std::endl;
	// std::cout << "\tvel->phi: " << vel->phi.ToString(10, 10) << std::endl;
	// std::cout << "\tequ: " << (-(gammaRThetaTheta * vel->theta * vel->theta) - (gammaRPhiPhi * vel->phi * vel->phi)).ToString(10, 10) << std::endl;

	floating gammaThetaPhiPhi	= -sin(pos->theta)*cos(pos->theta);
	floating gammaThetaRTheta	= 1/pos->r;
	floating gammaThetaThetaR	= floating(gammaThetaRTheta);

	floating gammaPhiRPhi		= floating(gammaThetaRTheta);
	floating gammaPhiPhiR		= floating(gammaThetaRTheta);
	floating gammaPhiThetaPhi	= cot(pos->theta);
	floating gammaPhiPhiTheta	= floating(gammaPhiThetaPhi);

	if (debug2) std::cout << " - christoffel symbols";

	// if (mpfr_get_exp(gammaTTR.num) > 0) std::cout << "gammaTTR: " << gammaTTR.ToString(10, 5) << std::endl;
	// if (mpfr_get_exp(gammaTRT.num) > 0) std::cout << "gammaTRT: " << gammaTRT.ToString(10, 5) << std::endl;
	// if (mpfr_get_exp(gammaRTT.num) > 0) std::cout << "gammaRTT: " << gammaRTT.ToString(10, 5) << std::endl;
	// if (mpfr_get_exp(gammaRRR.num) > 0) std::cout << "gammaRRR: " << gammaRRR.ToString(10, 5) << std::endl;
	/*if (mpfr_get_exp(gammaRThetaTheta.num) > 0)*/ std::cout << "gammaRThetaTheta: " << gammaRThetaTheta.ToString(10, 5) << std::endl;
	/*if (mpfr_get_exp(gammaRPhiPhi.num) > 0)*/ std::cout << "gammaRPhiPhi: " << gammaRPhiPhi.ToString(10, 5) << std::endl;
	// if (mpfr_get_exp(gammaThetaRTheta.num) > 0) std::cout << "gammaThetaRTheta: " << gammaThetaRTheta.ToString(10, 5) << std::endl;
	// if (mpfr_get_exp(gammaThetaThetaR.num) > 0) std::cout << "gammaThetaThetaR: " << gammaThetaThetaR.ToString(10, 5) << std::endl;
	// if (mpfr_get_exp(gammaThetaPhiPhi.num) > 0) std::cout << "gammaThetaPhiPhi: " << gammaThetaPhiPhi.ToString(10, 5) << std::endl;
	// if (mpfr_get_exp(gammaPhiRPhi.num) > 0) std::cout << "gammaPhiRPhi: " << gammaPhiRPhi.ToString(10, 5) << std::endl;
	// if (mpfr_get_exp(gammaPhiPhiR.num) > 0) std::cout << "gammaPhiPhiR: " << gammaPhiPhiR.ToString(10, 5) << std::endl;
	// if (mpfr_get_exp(gammaPhiThetaPhi.num) > 0) std::cout << "gammaPhiThetaPhi: " << gammaPhiThetaPhi.ToString(10, 5) << std::endl;
	// if (mpfr_get_exp(gammaPhiPhiTheta.num) > 0) std::cout << "gammaPhiPhiTheta: " << gammaPhiPhiTheta.ToString(10, 5) << std::endl;

	acceleration->t		= -(gammaTTR * vel->t * vel->r) - (gammaTRT * vel->r * vel->t);
	acceleration->r		= -(gammaRTT * vel->t * vel->t) - (gammaRRR * vel->r * vel->r) - (gammaRThetaTheta * vel->theta * vel->theta) - (gammaRPhiPhi * vel->phi * vel->phi);
	acceleration->theta	= -(gammaThetaRTheta * vel->r * vel->theta) - (gammaThetaThetaR * vel->theta * vel->r) - (gammaThetaPhiPhi * vel->phi * vel->phi);
	acceleration->phi	= -(gammaPhiRPhi * vel->r * vel->phi) - (gammaPhiPhiR * vel->phi * vel->r) - (gammaPhiThetaPhi * vel->theta * vel->phi) - (gammaPhiPhiTheta * vel->phi * vel->theta);

	if (debug2) std::cout << " - acceleration";
	// Deallocate cristoffel symbols
	gammaTTR.clear();
	gammaTRT.clear();
	gammaRTT.clear();
	gammaRRR.clear();
	gammaRThetaTheta.clear();
	gammaRPhiPhi.clear();
	gammaThetaRTheta.clear();
	gammaThetaThetaR.clear();
	gammaThetaPhiPhi.clear();
	gammaPhiRPhi.clear();
	gammaPhiPhiR.clear();
	gammaPhiThetaPhi.clear();
	gammaPhiPhiTheta.clear();
	if (debug2) std::cout << " - deallocated christoffel symbols";

	/*if (debug)*/ std::cout << "acceleration: " << vec3(*acceleration).toCartesian().ToString(10, 10) << std::endl;
	(*vel) += (*acceleration)*EPSILON;
	if (debug) std::cout << "vel: " << vec3(*vel).toCartesian().ToString(10, 10) << std::endl;
	(*pos) += (*vel);
	if (debug) std::cout << "pos: " << vec3(*pos).toCartesian().ToString(10, 10) << std::endl;
	if (debug2) std::cout << " - new position";
	// Write output to file instead of console
	vec3 vec = vec3(*pos);
	if (debug2) std::cout << " - created vec";
	vec = vec.toCartesian();
	if (debug2) std::cout << " - toCartesian";
	*outputFile << vec.ToString(10, 10) << std::endl;
	if (debug2) std::cout << " - wrote to file" << std::endl;
}

static void rk4(vec4* pos, vec4* vel, vec4* acceleration, floating mass, ofstream* outputFile, bool debug, bool debug2) {
	// k1 
}

int main() {
	vec4 pos = vec4(0.0, (vec3(0, 0, -300).toSphere()));
	pos.setName("pos");
	std::cout << "pos orig: " << vec3(pos).ToString(10, 10, false) << std::endl;
	vec4 vel = vec4(EPSILON, vec3(0,0,EPSILON).toSphere());
	vel.setName("vel");
	vec4 acceleration = vec4();
	acceleration.setName("acceleration");

	floating mass = floating(5);

	// Create and start writing to file
	ofstream outputFile("D:\\school\\SpaceScience\\blackHoleSim\\C_GNU\\3DOutput.txt");

	// The Schwarzschild Radius
	floating S_r = 2 * G * mass/(c*c);
	
	int iterations = 1000;
	for (int i = 0; i < iterations; i++) {
		step(&pos, &vel, &acceleration, mass, &outputFile, false, false);
		/*if (i % (iterations/10) == 0)*/ std::cout << i << std::endl;
		if (pos.r <= S_r) break;
	}
	


	// Close file stream
	outputFile.close();
	// Exit
    return 0;
}