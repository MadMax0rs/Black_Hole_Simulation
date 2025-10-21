
fov = 81
OriginalDistance = 100
camPos = vec3(100, 1.57, 0.57)//getPolar(vec3(0, -OriginalDistance, 0))
camDir = vec3(-0.999861, -0.0156978, -0.00569921)//normalize(-camPos)
camRight = normalize(cross(camDir, camUp))
ndc.x = 0//(400 - 800 / 2.0) / 800 * aspectRatio
ndc.y = 0//(fragCoord.y - iResolution.y / 2.0) / iResolution.y * aspectRatio
fovYRad = fovXRad = radians(fov)
camUpReal = cross(camDir, camRight)

normalize(
	camDir + (
		camRight * ndc.x * tan(fovXRad / 2.0)) + (
			camUpReal * ndc.y * tan(fovYRad / 2.0)
	)
) = 


(-0.743244, -0.668919, 0.0116689)