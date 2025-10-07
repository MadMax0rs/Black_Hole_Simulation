#version 330
precision highp float;
// https://chatgpt.com/s/t_68b3363462708191a45becc66da2e2ec

const vec2 size = vec2(5,5);
const float EPSILON = 0.001;
const float c = 3.0e8;
const float G = 6.674e-11;
const int maxSteps = 10000;
const float pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067;

uniform vec2 resolution;

out vec4 outColor;



vec3 CartesianToSpherical(vec3 cart) {
	// r, theta, phi
	// distance, azemuth (angle on the xy plane), polar angle(elevation)
	float theta = length(cart) == 0.0 ? 0.0 : acos(cart.z/length(cart));
	return vec3(
		length(cart),					// distance from 0,0,0
		theta,	// 0 -> north pole, pi/2 -> equator, pi -> south pole
		atan(cart.y, cart.x)			// 0 -> +x, pi/2 -> +z, -x -> pi, -z -> 3pi/2
	);
}
vec3 SphericalToCartesian(vec3 sphere) {
	// x, y, z
	// width (left-right), height, depth
	return vec3(
		sphere.x*sin(sphere.y) * cos(sphere.z),
		sphere.x*sin(sphere.y) * sin(sphere.z),
		sphere.x*cos(sphere.y)
	);
}

struct Ray {
	vec3 pos;
	vec3 dir;
	vec3 col;
	float b;
	vec2 momentum;
	int numSteps;

} ray;
struct RK4Components {
	// RK4 computents
	float t;
	float r;
	float theta;
	float phi;

};
RK4Components rk4Vector, rk4PVector;

RK4Components mul(RK4Components rk4, float a) {
	return RK4Components(
		rk4.t*a,
		rk4.r*a,
		rk4.theta*a,
		rk4.phi*a
	);
}
RK4Components mul(RK4Components rk4, RK4Components rk4_2) {
	return RK4Components(
		rk4.t*rk4_2.t,
		rk4.r*rk4_2.r,
		rk4.theta*rk4_2.theta,
		rk4.phi*rk4_2.phi
	);
}
RK4Components div(RK4Components rk4, float a) {
	return RK4Components(
		rk4.t/a,
		rk4.r/a,
		rk4.theta/a,
		rk4.phi/a
	);
}
RK4Components div(RK4Components rk4, RK4Components rk4_2) {
	return RK4Components(
		rk4.t/rk4_2.t,
		rk4.r/rk4_2.r,
		rk4.theta/rk4_2.theta,
		rk4.phi/rk4_2.phi
	);
}
RK4Components add(RK4Components rk4, float a) {
	return RK4Components(
		rk4.t+a,
		rk4.r+a,
		rk4.theta+a,
		rk4.phi+a
	);
}
RK4Components add(RK4Components rk4, RK4Components rk4_2) {
	return RK4Components(
		rk4.t+rk4_2.t,
		rk4.r+rk4_2.r,
		rk4.theta+rk4_2.theta,
		rk4.phi+rk4_2.phi
	);
}
RK4Components add(RK4Components rk4_1, RK4Components rk4_2, RK4Components rk4_3, RK4Components rk4_4) {
	return add(add(add(rk4_1, rk4_2), rk4_3), rk4_4);
}
float cot(float theta) {
	return cos(theta)/sin(theta);
}

vec3 point = vec3(0,0,0);
float radius = 1.0;
float mass = 6.74258316e26;		// (c*c*radius)/(2G)


RK4Components[2] RK4(RK4Components rk4, RK4Components p) {
	// f(r)
	float fOfR = 1.0 - radius/rk4.r;
	// f'(r)
	float fPrime = radius/(rk4.r*rk4.r);

	float gammaTTR, gammaTRT = fPrime / (2.0 * fOfR);
	float gammaRTT = fOfR * fPrime * c*c / 2.0;
	float gammaRRR = -2.0 * fOfR / fPrime;
	float gammaRThetaTheta = -rk4.r * fOfR;
	float gammaRPhiPhi = -rk4.r * fOfR * sin(rk4.theta)*sin(rk4.theta);
	float gammaThetaRTheta, gammaThetaThetaR = 1.0/rk4.r;
	float gammaThetaPhiPhi = -sin(rk4.theta)*cos(rk4.theta);
	float gammaPhiRPhi, gammaPhiPhiR = 1.0/rk4.r;
	float gammaPhiThetaPhi, gammaPhiPhiTheta = cot(rk4.theta);

	float F_1 = p.t;
	float F_2 = p.r;
	float F_3 = p.theta;
	float F_4 = p.phi;
	float F_5 = -2.0 * gammaTTR * p.t * p.r;
	float F_6 =
		-gammaRTT * p.t*p.t +
		-gammaRRR * p.r*p.r +
		-gammaRThetaTheta * p.theta*p.theta +
		-gammaRPhiPhi * p.phi*p.phi;
	float F_7 = -2.0 * gammaThetaRTheta * p.r * p.theta
		- gammaThetaPhiPhi * p.phi * p.phi;
	float F_8 = -2.0 * gammaPhiRPhi * p.r * p.theta
		- 2.0 * gammaPhiThetaPhi * p.theta * p.phi;

	RK4Components[2] ret;
	ret[0] = RK4Components(F_1, F_2, F_3, F_4);
	ret[1] = RK4Components(F_5, F_6, F_7, F_8);
	return ret;
}

bool stepRay() {
	ray.numSteps++;
	

	RK4Components[2] k_1 = RK4(rk4Vector, rk4PVector);
	RK4Components[2] k_2 = RK4(
		add(rk4Vector, mul(k_1[0], EPSILON/2.0)),
		add(rk4Vector, mul(k_1[1], EPSILON/2.0))
	);
	RK4Components[2] k_3 = RK4(
		add(rk4Vector, mul(k_2[0], EPSILON/2.0)),
		add(rk4Vector, mul(k_2[1], EPSILON/2.0))
	);
	RK4Components[2] k_4 = RK4(
		add(rk4Vector, mul(k_3[0], EPSILON/2.0)),
		add(rk4Vector, mul(k_3[1], EPSILON/2.0))
	);

	//ray = add(ray, mul(add(k_1, mul(k_2, 2.0), mul(k_3, 2.0), k_4), (EPSILON/6.0)));


	// Step ray
	//ray.pos += ray.dir*EPSILON;

	// If ray hits something, return true
	if (rk4Vector.r < radius) {
		// Hit event horizon
		ray.col = vec3(0);
		return true;
	} else if (rk4Vector.r > 10.0) {
		// Go into space
		ray.col = vec3(0.7);
		return true;
	}
	if (ray.numSteps > maxSteps) {
		ray.col = vec3(1, 0, 0);
		return true;
	}

	// if hits sphere behind black hole (coordinates relative to black hole)
	if (length(SphericalToCartesian(vec3(rk4Vector.r, rk4Vector.theta, rk4Vector.phi)) - vec3(0,0,2)) < 0.1) {
		ray.col = vec3(0,0,1);
		return true;
	}
	return false;
}

bvec3 or(bvec3 a, bvec3 b) {
	return bvec3(a.x || b.x, a.y || b.y, a.z || b.z);
}


vec2 worldToScreenCenter(vec3 cartesianWorld) {
	float aspectRatio = resolution.x/resolution.y;
	vec2 normUV = gl_FragCoord.xy/resolution*vec2(aspectRatio, 1);

	return cartesianWorld.zx/vec2(10.0/aspectRatio,10.0) - normUV + vec2(aspectRatio/2.0, 0.5);
}

void main() {
	// Init Ray
	//ray.pos = vec3((gl_FragCoord.xy-resolution/2.0)/(resolution.x/3.0), 0);
	ray.pos = vec3(0,0,-5);
	ray.dir = vec3(0, 0, 1);
	ray.col = vec3(1);
	// Calculate impact parameter (𝑏)
	ray.b = length(cross(ray.dir, (ray.pos - point)))/length(ray.dir);

	vec3 coords = CartesianToSpherical(ray.pos - point);
	rk4Vector = RK4Components(0.0, coords.x, coords.y, coords.z);


	vec3 newCoords = (CartesianToSpherical(point - (ray.pos + vec3(0,0,EPSILON))) - coords)/EPSILON;

	float fOfr_0 = 1.0 - radius/coords.x;
	float rSquare = coords.x*coords.x;
	float pt = sqrt(
		// (f*(p^r)^2+
		(newCoords.x*newCoords.x/fOfr_0)+
		// (r_0)^2(p^θ)^2
		(rSquare*newCoords.y*newCoords.y)+
		// (r_0)^2*sin(θ)^2(p^ϕ)^2)/
		(rSquare * sin(coords.y)*sin(coords.y) * newCoords.z*newCoords.z) /
		// f*c^2
		fOfr_0*c*c
	);

	rk4PVector = RK4Components(pt, coords.x, coords.y, coords.z);

	bool hit = false;

	float aspectRatio = resolution.x/resolution.y;
	vec2 normUV = gl_FragCoord.xy/resolution*vec2(aspectRatio, 1);

	outColor = vec4(normUV, 0, 1);
	// return;
	while (!hit) {
		// Trace Ray
		vec3 vector = SphericalToCartesian(vec3(rk4Vector.r,rk4Vector.theta,rk4Vector.phi));

		// Check vector is valid
		// bvec3 invalid = or(isnan(vector), isinf(vector));
		// if (invalid.x || invalid.y || invalid.z) {
		// 	outColor = vec4(1,0,1,1);
		// 	return;
		// }

		// Trace ray path
		if (length(worldToScreenCenter(vector)) < 0.005) {
			vec3 numSteps = vec3(float(ray.numSteps)/1.0);
			outColor = vec4(0,0,1,1);
			return;
		}

		hit = stepRay();
	}

	vec3 numSteps = vec3(float(ray.numSteps)/15000.0);
	//outColor = vec4(1.0);
}