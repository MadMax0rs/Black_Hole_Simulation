import math
import numpy as np

class vec3:
	x: float
	y: float
	z: float

	def __init__(self, x: float = 0, y: float = 0, z: float = 0):
		self.x = x
		self.y = y
		self.z = z

	def __str__(self):
		return f"({self.x},{self.y},{self.z})"

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
	
	def __sub__(self, other):
		if (type(other) != vec3):
			raise TypeError
		return vec3(self.x - other.x, self.y - other.y, self.z - other.z)

def CartesianToSphereical(cart: vec3) -> vec3:
	# r, theta, phi
	# distance, azemuth (angle on the xy plane), zenith (elevation)
	theta: float = 0.0
	try:
		theta = math.acos(cart.z/cart.length)
	except ZeroDivisionError:
		theta = 0
		pass

	return vec3(
		cart.length,					# distance from 0,0,0
		theta,	# 0 -> north pole, pi/2 -> equator, pi -> south pole
		math.atan2(cart.y, cart.x)		# 0 -> +x, pi/2 -> +z, -x -> pi, -z -> 3pi/2
	)

def SphericalToCartesian(sphere: vec3) -> vec3:
	# x, y, z
	# width (left-right), height, depth
	return vec3(
		sphere.r*math.sin(sphere.theta) * math.cos(sphere.phi),
		sphere.r*math.sin(sphere.theta) * math.sin(sphere.phi),
		sphere.r*math.cos(sphere.theta)
	)

def test(cartesian: vec3):
	newCoords: vec3 = CartesianToSphereical(cartesian)
	newCoords = SphericalToCartesian(newCoords)
	print("Error:", (newCoords-cartesian).length)
	return


test(vec3(0,0,0))
# test(vec3(2,2,2))
# test(vec3(1,2,3))
# test(vec3(math.pi,0.0001,1000000))
# test(vec3(9467e10,999999999,1000000))
# test(vec3(-15,12,-12345))