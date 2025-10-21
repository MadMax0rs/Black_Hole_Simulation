#version 330

uniform vec2 uResolution;
uniform float uTime;

out vec4 fragColor;

#define OriginalDistance 100.0
float M = 5.0;

// Fork of "Black Hole RayTracer/RayMarcher" by ShaderGott420. https://shadertoy.com/view/4XsGWj
// 2024-02-15 13:38:09


float PI = 3.14159265358979323846;
float diskThickness = 0.01;								// Thickness of the disk in radians
float equatorialPlane = 3.14159265358979323846 / 2.0;	// Equatorial plane angle in radians


vec3 getCartesian(vec3 polar)
{
	float r = polar.x;
	float theta = polar.y;
	float phi = polar.z;

	return vec3(
		r * sin(theta) * cos(phi),
		r * cos(theta),
		r * sin(theta) * sin(phi)
	);
}
vec3 getPolar(vec3 cartesian)
{
	float r = length(cartesian);				// GLSL's length function computes the magnitude of the vector
	float theta = acos(cartesian.z / r);		// GLSL uses .x, .y, .z for vector components
	float phi = atan(cartesian.y, cartesian.x);	// atan(y, x) is used instead of atan2

	return vec3(r, theta, phi);
}
// Stars texture from Kali
// https://www.shadertoy.com/view/XlfGRj
// Very Pretty

void assignChristoffelSymbols(vec3 polarCoords, out float trt, out float ttr, out float rrr, out float rtt, out float rphiphi, out float rthetatheta, out float thetartheta, out float thetathetar, out float thetaphiphi, out float phirphi, out float phiphir, out float phithetaphi, out float phiphitheta)
{
	float r = polarCoords.x;     // Radial distance
	float theta = polarCoords.y; // Polar angle
	float phi = polarCoords.z;   // Azimuthal angle (unused in the current symbols but included for completeness)

	// Now assign each Christoffel symbol based on r, theta, and predefined M
	trt = M / (r * (r - 2.0 * M));
	ttr = M / (r * (r - 2.0 * M));
	rrr = -M / (r * (r - 2.0 * M));
	rtt = (M * (r - 2.0 * M)) / pow(r, 3.0);
	rphiphi = -1.0 * (r - 2.0 * M) * pow(sin(theta), 2.0);
	rthetatheta = -1.0 * (r - 2.0 * M);
	thetartheta = 1.0 / r;
	thetathetar = 1.0 / r;
	thetaphiphi = -sin(theta) * cos(theta);
	phirphi = 1.0 / r;
	phiphir = 1.0 / r;
	phithetaphi = 1.0 / tan(theta);
	phiphitheta = 1.0 / tan(theta);
}
vec3 updateVelocity(vec3 VE, vec3 P)
{
	float r = P.x;
	float theta = P.y;
	float phi = P.z;

	// Variables to hold the Christoffel symbols
	float trt, ttr, rrr, rtt, rphiphi, rthetatheta, thetartheta, thetathetar, thetaphiphi, phirphi, phiphir, phithetaphi, phiphitheta;

	// Assuming assignChristoffelSymbols is defined to calculate and assign Christoffel symbols
	assignChristoffelSymbols(P, trt, ttr, rrr, rtt, rphiphi, rthetatheta, thetartheta, thetathetar, thetaphiphi, phirphi, phiphir, phithetaphi, phiphitheta);

	float vr = VE.x;
	float vtheta = VE.y;
	float vphi = VE.z;

	float dv_r = -1.0 * rrr * (vr * vr) - rphiphi * (vphi * vphi) - rthetatheta * (vtheta * vtheta);
	float dv_theta = -1.0 * thetartheta * (vr * vtheta) - thetathetar * (vtheta * vr) - thetaphiphi * (vphi * vphi);
	float dv_phi = -1.0 * phirphi * (vr * vphi) - phiphir * (vphi * vr) - phithetaphi * (vtheta * vphi) - phiphitheta * (vphi * vtheta);

	vec3 dVE = vec3(dv_r, dv_theta, dv_phi);
	VE += dVE;

	return VE;
}

vec3 sphereCenter = vec3(0.0, 0.0, 0); // Sphere at origin

void main()
{
    vec2 uv = gl_FragCoord.xy/uResolution.xy;

	float sphereRadius = 2.0 * M * 16.0;

	float aspectRatio = uResolution.x / uResolution.y;

	vec3 camPos = getPolar(vec3(0, -OriginalDistance, 0));//vec3((OriginalDistance - (Speed * iTime)) * sin(iTime * camSpeed), 0, (OriginalDistance - (Speed * iTime)) * cos(iTime * camSpeed)); // Position of the camera


	// Calculate camera basis vectors inside mainImage
	//vec3 camDir = normalize(-camPos);

	// Calculate the ray direction from camera parameters
	// normalize(normalize(-getPolar(vec3(0, -OriginalDistance, 0))) + (normalize(cross(normalize(-getPolar(vec3(0, -OriginalDistance, 0))), vec3(0.0, 0.0, 1.0))) * ((fragCoord.x - iResolution.x / 2.0) / iResolution.x * aspectRatio) * tan(radians(fov) / 2.0)) + (cross(normalize(-getPolar(vec3(0, -OriginalDistance, 0))), normalize(cross(normalize(-getPolar(vec3(0, -OriginalDistance, 0))), vec3(0.0, 0.0, 1.0)))) * ((fragCoord.y - iResolution.y / 2.0) / iResolution.y) * tan(radians(fov) / 2.0)))
	vec2 normUV = uv - 0.5;
	normUV.x *= aspectRatio;
	vec2 temp = vec2(normUV.x, -normUV.y);
	vec3 rayDir = vec3(-1, 0, 0);//vec3(-sqrt(1-temp.x*temp.x-temp.y*temp.y), temp);//normalize(camDir + (camRight * ndc.x * tan(fovXRad / 2.0)) + (camUpReal * ndc.y * tan(fovYRad / 2.0)));

	// Use this ray direction with the camera's position to cast rays into the scene
	vec3 rayOrigin = getPolar(vec3(-OriginalDistance, normUV.x*1050.0, -normUV.y*1050.0));//camPos;

	// Use this ray direction with the camera's position to cast rays into the scene
	vec3 rayMarchedPosition = camPos;
	float marchDistance = 0.0;
	// Sphere definition
	float dist = length(rayOrigin);
	float Precision = smoothstep(20.0 * M, 0.0, dist) * 200.0 + (smoothstep(OriginalDistance, 20.0 * M, dist) + 1.0) * 5.0 * M;
	
	vec3 rayMarchedVelocity = getPolar(rayMarchedPosition + rayDir / 20.0) - getPolar(rayMarchedPosition); // the divisor of rayDir is dependent on a few things




	// MAIN LOOP
	
	bool noHit = true;
	float OuterDiskRadius = 6.0 * M;

	// Prolly gon do some polar Transformations.
	rayMarchedPosition = getPolar(rayMarchedPosition);
	for (int i = 2; i < int(900.0 * Precision); i++)
	{
		// Hit Event Horizon
		if (rayMarchedPosition.x <= 2.1 * M)
		{
			noHit = false;
			fragColor = vec4(0, 0, 0, 1);
			break;
		}
		
		// Shorten render time by having a cutoff for gravity's affect
		/*if (dot(getCartesian(rayMarchedVelocity), -getCartesian(rayMarchedPosition)) < 0.0 && rayMarchedPosition.x > effectiveRadius)
		{
			break;
		}*/

		// Hit accretian disk
		if (rayMarchedPosition.x <= OuterDiskRadius && abs(rayMarchedPosition.y - equatorialPlane) <= diskThickness)
		{
			vec4 text = vec4(1);//texture(iChannel0, RMP.zx);
			if (dot(normalize(text), vec4(0.000, 0.000, 0.000, 1.0)) < 0.95)
			{
				fragColor = text;
				return;
			}
		};

		// March
		rayMarchedVelocity = updateVelocity(rayMarchedVelocity, rayMarchedPosition);
		rayMarchedPosition += rayMarchedVelocity;
	}
	// Didn't hit
	if (noHit)
	{
        fragColor = vec4(uv, 0.0, 1.0);
	}
}
