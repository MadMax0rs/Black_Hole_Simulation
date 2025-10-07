from math import *
import numpy as np
import pygame as pg
import sys

EPSILON: float = 0.001
c: float = 1#3.0e8
G: float = 1#6.674e-11
MAX_STEPS: float = 10000


print("tan(pi) = ", tan(pi))

class vec3:
	x: float
	y: float
	z: float

	def __init__(self, x: float = 0, y: float = 0, z: float = 0):
		self.x = x
		self.y = y
		self.z = z

	def __str__(self):
		return f"({self.x},\n{self.y},\n{self.z})"
	
	def __add__(self, other):
		if (type(other) != vec3):
			raise TypeError
		return vec3(self.x + other.x, self.y + other.y, self.z + other.z)
	def __sub__(self, other):
		if (type(other) != vec3):
			raise TypeError
		return vec3(self.x - other.x, self.y - other.y, self.z - other.z)
	def __truediv__(self, other):
		if (type(other) != float):
			raise TypeError
		return vec3(self.x / other, self.y / other, self.z / other)

	# Multiplication: vec3 * scalar
	def __mul__(self, other):
		if isinstance(other, (int, float)):
			return vec3(
				self.x * other,
				self.y * other,
				self.z * other
			)
		if isinstance(other, (vec3)):
			return vec3(
				self.x * other.x,
				self.y * other.y,
				self.z * other.z
			)
		raise TypeError("Can only multiply vec3 by a scalar")

	# Right-side multiplication: scalar * vec3
	def __rmul__(self, other):
		return self.__mul__(other)

	def __repr__(self):
		return f"Vec3(x={self.x}, y={self.y}, z={self.z})"

	@property
	def r(self):
		return self.x
	@property
	def theta(self):
		return self.y
	@property
	def phi(self):
		return self.z
	
	@property
	def length(self):
		return sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
	


class vec4:
	def __init__(self, *args):
		if len(args) == 4 and all(isinstance(a, (int, float)) for a in args):
			# Constructor: Vec4(t, r, theta, phi)
			self.t, self.r, self.theta, self.phi = args
		elif len(args) == 2 and isinstance(args[0], (int, float)) and isinstance(args[1], vec3):
			# Constructor: Vec4(t, Vec3)  (map x→r, y→theta, z→phi)
			t, v3 = args
			self.t = t
			self.r, self.theta, self.phi = v3.x, v3.y, v3.z
		else:
			raise TypeError("Vec4 must be initialized with (t, r, theta, phi) or (t, Vec3)")
	def __str__(self) -> str:
		return f"(t={self.t}, r={self.r}, θ={self.theta}, φ={self.phi})"
	
	def __repr__(self) -> str:
		return f"Vec4(t={self.t}, r={self.r}, θ={self.theta}, φ={self.phi})"

	# Addition: vec4 + vec4
	def __add__(self, other):
		if isinstance(other, vec4):
			return vec4(
				self.t + other.t,
				self.r + other.r,
				self.theta + other.theta,
				self.phi + other.phi
			)
		raise TypeError("Can only add vec4 to vec4")

	# Subtraction: vec4 - vec4
	def __sub__(self, other):
		if isinstance(other, vec4):
			return vec4(
				self.t - other.t,
				self.r - other.r,
				self.theta - other.theta,
				self.phi - other.phi
			)
		raise TypeError("Can only subtract vec4 from vec4")

	# Multiplication: vec4 * scalar
	def __mul__(self, other):
		if isinstance(other, (int, float)):
			return vec4(
				self.t * other,
				self.r * other,
				self.theta * other,
				self.phi * other
			)
		raise TypeError("Can only multiply vec4 by a scalar")

	# Right-side multiplication: scalar * vec4
	def __rmul__(self, other):
		return self.__mul__(other)

	# True division: vec4 / scalar
	def __truediv__(self, other):
		if isinstance(other, (int, float)):
			return vec4(
				self.t / other,
				self.r / other,
				self.theta / other,
				self.phi / other
			)
		raise TypeError("Can only divide vec4 by a scalar")
	
	@property
	def yzw(self) -> vec3:
		return vec3(self.r, self.theta, self.phi)
	

class State:
	def __init__(self, rk4: vec4, p: vec4):
		self.rk4: vec4 = rk4
		self.p: vec4 = p

	# Multiplication: State * scalar
	def __mul__(self, other):
		if isinstance(other, (int, float)):
			return State(
				self.rk4 * other,
				self.p * other,
			)
		raise TypeError("Can only multiply vec4 by a scalar")

	# Right-side multiplication: scalar * vec4
	def __rmul__(self, other):
		return self.__mul__(other)
	
	# Addition: State + State
	def __add__(self, other):
		if isinstance(other, State):
			return State(
				self.rk4 + other.rk4,
				self.p + other.p,
			)
		raise TypeError("Can only add State to State")
	
	def split(self) -> tuple[vec4, vec4]:
		return (self.rk4, self.p)
	
class Ray:
	def __init__(self, pos: vec3, dir: vec3, col: vec3, numSteps: float = 0):
		self.pos: vec3 = pos
		self.dir: vec3 = dir
		self.col: vec3 = col
		self.numSteps: float = numSteps

def CartesianToSpherical(cart: vec3) -> vec3:
	# r, theta, phi
	# distance, azemuth (angle on the xy plane), zenith (elevation)
	theta: float = 0.0
	try:
		theta = acos(cart.z/cart.length)
	except ZeroDivisionError:
		theta = 0
		pass

	return vec3(
		cart.length,					# distance from 0,0,0
		theta,	# 0 -> north pole, pi/2 -> equator, pi -> south pole
		atan2(cart.y, cart.x)		# 0 -> +x, pi/2 -> +z, -x -> pi, -z -> 3pi/2
	)

def SphericalToCartesian(sphere: vec3) -> vec3:
	# x, y, z
	# width (left-right), height, depth
	return vec3(
		sphere.r*sin(sphere.theta) * cos(sphere.phi),
		sphere.r*sin(sphere.theta) * sin(sphere.phi),
		sphere.r*cos(sphere.theta)
	)

def cot(theta: float) -> float:
	theta = max(EPSILON, min(pi - EPSILON, theta))
	return 1.0 / tan(theta)

def RK4(state: State) -> State:
	rk4, p = state.split()
	# f(r)
	fOfR: float = 1.0 - radius/rk4.r
	# f'(r)
	fPrime: float = radius/(rk4.r*rk4.r)
	test = rk4.r * rk4.r
	gammaTTR = gammaTRT = fPrime / (2.0 * fOfR)
	gammaRTT: float = fOfR * fPrime * c*c / 2.0
	gammaRRR: float = fPrime / (-2.0 * fOfR)
	gammaRThetaTheta: float = -rk4.r * fOfR
	gammaRPhiPhi: float = -rk4.r * fOfR * sin(rk4.theta)*sin(rk4.theta)
	gammaThetaRTheta = gammaThetaThetaR = 1.0/rk4.r
	gammaThetaPhiPhi: float = -sin(rk4.theta)*cos(rk4.theta)
	gammaPhiRPhi = gammaPhiPhiR = 1.0/rk4.r
	gammaPhiThetaPhi = PhiPhiTheta = cot(rk4.theta)

	# (t, r, theta, phi), (-2gammaTTR*t*r, -gammaRTT*t^2+-gammaRRR*r^2+-gammaRThetaTheta)
	F_1: float = p.t
	F_2: float = p.r
	F_3: float = p.theta
	F_4: float = p.phi
	F_5: float = -2.0 * gammaTTR * p.t * p.r
	F_6: float = (
		-gammaRTT * p.t*p.t +
		-gammaRRR * p.r*p.r +
		-gammaRThetaTheta * p.theta*p.theta +
		-gammaRPhiPhi * p.phi*p.phi
	)
	F_7: float = (
		-2.0 * gammaThetaRTheta * p.r * p.theta
		- gammaThetaPhiPhi * p.phi * p.phi
	)
	F_8: float = (
		-2.0 * gammaPhiRPhi * p.r * p.theta
		- 2.0 * gammaPhiThetaPhi * p.theta * p.phi
	)

	return State(vec4(F_1, F_2, F_3, F_4), vec4(F_5, F_6, F_7, F_8))


point: vec3 = vec3(0,0,0)
radius: float = 1
mass: float = 6.74258316e26

camPos: vec3 = vec3(0,0,-5)
ray: Ray = Ray(vec3(0,0,-5), vec3(0,0,1), vec3())

coords: vec3 = CartesianToSpherical(ray.pos - point)
rk4Vector: vec4 = vec4(0, coords)


# Calculate inputs for the p vector
Sphericalvelocity: vec3 = (CartesianToSpherical((ray.dir - point) * EPSILON))

fOfr_0: float = 1.0 - radius/coords.x
rSquare: float = coords.x**2
pt: float = sqrt(
	# (f*(p^r)^2+
	(Sphericalvelocity.r**2/fOfr_0 +
	# (r_0)^2(p^θ)^2
	rSquare*Sphericalvelocity.theta**2 +
	# (r_0)^2*sin(θ)^2(p^ϕ)^2)/
	rSquare * sin(coords.theta)**2 * Sphericalvelocity.phi**2) /
	# f*c^2
	fOfr_0*c**2
)

rk4PVector = vec4(pt, coords)



def stepRay(rk4Vector: vec4, rk4PVector: vec4) -> tuple[bool, State]:
	ray.numSteps += 1
	
	state: State = State(rk4Vector, rk4PVector)
	k_1: State = RK4(state)
	k_2: State = RK4(state + k_1 * (EPSILON/2.0))
	k_3: State = RK4(state + k_2 * (EPSILON/2.0))
	k_4: State = RK4(state + k_3 * EPSILON)

	state: State = State(rk4Vector, rk4PVector)
	state += (k_1 + 2*k_2 + 2*k_3 + k_4) * (EPSILON/6.0)



	# Step ray
	#ray.pos += ray.dir*EPSILON

	# If ray hits something, return true
	if (rk4Vector.r < radius):
		# Hit event horizon
		ray.col = vec3(0)
		return (True, state)
	elif (rk4Vector.r > 10.0):
		# Go into space
		ray.col = vec3(0.7)
		return (True, state)
	if (ray.numSteps > MAX_STEPS):
		ray.col = vec3(1, 0, 0)
		return (True, state)

	# if hits sphere behind black hole (coordinates relative to black hole)
	if (SphericalToCartesian(vec3(rk4Vector.r, rk4Vector.theta, rk4Vector.phi) - vec3(0,0,2)).length < 0.1):
		ray.col = vec3(0,0,1)
		return (True, state)
	
	return (False, state)

def WorldToScreen(world: vec3) -> tuple[float, float]:
	pass

pg.init()

screen_width = 800
screen_height = 600
screen: pg.Surface = pg.display.set_mode((screen_width, screen_height))
stop = False
while not stop:
	for event in pg.event.get():
		if (event.type == pg.QUIT):
			stop = True

	running, state = stepRay(rk4Vector, rk4PVector)
	rk4Vector, rk4PVector = state.split()
	print(rk4Vector.yzw)
	pg.draw.circle(screen, (1,0,0), (SphericalToCartesian()))

pg.quit()