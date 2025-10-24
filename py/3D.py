import math

EPSILON: float = 0.00001
def strFunc(string: str) -> str:
	if (string.__contains__("e")):
		return string.replace("e", "*10^{") + "}"
	return string
class vec3:
	x: float
	y: float
	z: float

	def __init__(self, x: float = 0, y: float = 0, z: float = 0):
		self.x = x
		self.y = y
		self.z = z

	def toSphere(self):
		if self.x == 0 and self.y == 0 and self.z == 0:
			return vec3(0,0,0)
		return vec3(
			self.length, 
			math.acos(self.y/self.length), 
			math.atan2(self.z, self.x)
		)

	def toCart(self):
		return vec3(
			self.r * math.sin(self.theta) * math.cos(self.phi),
			self.r * math.cos(self.theta),
			self.r * math.sin(self.theta) * math.sin(self.phi)
		)

	def XZString(self):
		return f"({strFunc(str(self.x))}, {strFunc(str(self.z))})"
	def YZString(self):
		return f"({strFunc(str(self.y))}, {strFunc(str(self.z))})"
	def XYString(self):
		return f"({strFunc(str(self.x))}, {strFunc(str(self.y))})"
	def XZYString(self):

		return f"({strFunc(str(self.x))}, {strFunc(str(self.z))}, {strFunc(str(self.y))})"

	def __str__(self):
		return f"({strFunc(str(self.x))},{strFunc(str(self.y))},{strFunc(str(self.z))})"
	
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
		return math.sqrt(self.x*self.x+self.y*self.y+self.z*self.z)
	


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
	

saves = [
	{"pos": vec4(0, vec3(0, -2, -10).toSphere()), "vel": vec4(0, vec3(-0.01, 0.00024, 0)), "mass": 1, "rk4": False},
	{"pos": vec4(0, vec3(0, -2, -10).toSphere()), "vel": vec4(0, vec3(-0.01, 0.00023, 0.0001)), "mass": 1, "rk4": False},
	{"pos": vec4(0, vec3(5, -2, -10).toSphere()), "vel": vec4(0, vec3(-0.01, 0.01, 0.001)), "mass": 1, "rk4": True},
]
pos = vec4(0, vec3(5, -2, -10).toSphere())
vel = vec4(0, vec3(-0.01, 0.01, 0.001))
useRk4 = False
mass: float = 1

def load(index: int):
	global pos, vel, mass, useRk4
	pos = saves[index]["pos"]
	vel = saves[index]["vel"]
	mass = saves[index]["mass"]
	useRk4 = saves[index]["rk4"]

load(1)

def step(pos, vel) -> vec4:

	# Calculate christoffel symbols
	common = pos.r - 2.0*mass 

	ttr = trt = mass / (pos.r * common)
	
	# rtt: float = mass*common / (pos.r**3)
	rrr: float = -ttr
	rThetaTheta: float = -common
	rPhiPhi: float = rThetaTheta*(math.sin(pos.theta)**2)

	thetaPhiPhi: float = -math.sin(pos.theta)*math.cos(pos.theta)
	thetaRTheta = thetaThetaR = 1.0/pos.r

	phiRPhi = phiPhiR = thetaRTheta
	phiPhiTheta = phiThetaPhi = 1.0/math.tan(pos.theta)

	acc: vec4 = vec4(
		0,
		-rrr*(vel.r**2) - rThetaTheta*(vel.theta**2) - rPhiPhi*(vel.phi**2),
		-thetaPhiPhi*(vel.phi**2) - thetaRTheta*(vel.r*vel.theta) - thetaThetaR*(vel.theta*vel.r),
		-phiRPhi*(vel.r*vel.phi) - phiPhiR*(vel.phi*vel.r) - phiPhiTheta*(vel.phi*vel.theta) - phiThetaPhi*(vel.theta*vel.phi)
	)
	return acc

def rk4():
	global pos, vel
	p = pos
	v = vel
	k1 = step(p, v)
	v += k1
	p += v
	k2 = step(p * EPSILON/2, v * EPSILON/2)
	v += k1
	p += v
	k3 = step(p * EPSILON/2, v * EPSILON/2)
	v += k1
	p += v
	k4 = step(p * EPSILON, v * EPSILON)

	vel += (k1 + 2*k2 + 2*k3 + k4) * EPSILON/6
	pos += vel



ITERATIONS = 20000


with open("D:/school/SpaceScience/blackHoleSim/py/3DOutput.txt", "w") as file:

	for i in range(ITERATIONS):
		if (i%10 == 0):
			file.write((str(pos.yzw.toCart().XZYString()) + "\n"))
		if (i%1000 == 0):
			print(i)

		if (useRk4):
			rk4()
		else:
			acc = step(pos, vel)
			vel += acc
			pos += vel
		if pos.r < 2.0*mass:
			file.write(pos.yzw.toCart().XZString())
			break