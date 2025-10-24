#include <mpfr.h>
#include <vector>
#include <iostream>


PRECISION: int = 256
EPSILON: float = 0.0001
c: float = 3.0e8
G: float = 6.674e-11

class vec3:
	# (x,y,z), (t,r,phi)

	def __init__(self, x: float, y: float, z: float, precision: int = PRECISION):
		self.x = x
		self.y = y
		self.z = z

	def printToConsole(self):
		print(f"({self.x}, {self.y}, {self.z})")


	def __add__(self, other):
		return vec3(self.x + other.x, self.y + other.y, self.z + other.z)
	def __sub__(self, other):
		return vec3(self.x - other.x, self.y - other.y, self.z - other.z)
	def __mul__(self, other):
		return vec3(self.x * other.x, self.y * other.y, self.z * other.z)
	def __div__(self, other):
		return vec3(self.x / other.x, self.y / other.y, self.z / other.z)



def step(pos: vec3, vel: vec3, mass: float) -> vec3:

	
	# gammaRTT = GM/r^2(1-2GM/(rc^2))
	gammaRTT = G*mass/(pos.y*pos.y)


	# gammaTTR = GM/(r^2c^2)

	# gammaRRR = gammaTRT = gammaTTR = GM/(r^2c^2)(1/(1-2GM/(rc^2)))

	#gammaRRR = -GM/(r^2c^2)(1/(1-2GM/(rc^2)))
	mpfr_mul_d(gammaRRR, gammaRRR, -1, MPFR_RNDN)

	vec3 acceleration = vec3(0.0, 0.0, 0.0)
	# acceleration.x = (gammaTTR + gammaTRT)vel.t * vel.r
	mpfr_add(acceleration.x, gammaTTR, gammaTRT, MPFR_RNDN)
	mpfr_mul(acceleration.x, acceleration.x, vel.x, MPFR_RNDN)
	mpfr_mul(acceleration.x, acceleration.x, vel.y, MPFR_RNDN)

	# acceleration.y = gammaRTT * vel.t^2 - gammaRRR * vel.r^2
	mpfr_sqr(temp,vel.x , MPFR_RNDN)
	mpfr_mul(temp, temp, gammaRTT, MPFR_RNDN)
	mpfr_sqr(acceleration.y, vel.y, MPFR_RNDN)
	mpfr_mul(acceleration.y, acceleration.y, gammaRRR, MPFR_RNDN)
	mpfr_sub(acceleration.y, temp, acceleration.y, MPFR_RNDN)
	return acceleration
}

int main():


	vec3 pos = vec3(0, -20, 0)
	vec3 vel = vec3(0, c, 0)
	mpfr_t mass
	mpfr_set_d(mass, 1.0, MPFR_RNDN)

	vel += step(pos, vel, mass)
	pos += vel

	
	pos.print()

	# Exit
	return 0
}