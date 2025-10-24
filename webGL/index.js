const vertices = new Float32Array([
	-1,-1,
	-1, 1,
	1, 1,
	-1,-1,
	1, 1,
	1,-1
])

var program, gl

var canvas = document.createElement('canvas')
var ctx = canvas.getContext("2d")

window.onresize = reload

async function loadShaderFiles() {
	vertShaderSource = await (await fetch('http://localhost:8000/vert.glsl?v=' + Date.now())).text()
	fragShaderSource = await (await fetch('http://localhost:8000/frag.glsl?v=' + Date.now())).text()
	reload()
}

var vertShaderSource = ""
var fragShaderSource = ""

var canvas = document.createElement('canvas')
canvas.id = "canvas"
document.body.appendChild(canvas)

var ctx = canvas.getContext("2d")

loadShaderFiles()


function reload() {
	console.clear()
	canvas.remove()
	canvas = document.createElement('canvas')
	canvas.id = "canvas"
	// Debug window scale and aspect ratio
	console.log(`canvas ${canvas.width}x${canvas.height}	${canvas.width/canvas.height}`);					// Canvas is uneffected by dev. tools window
	console.log(`window ${window.innerWidth}x${window.innerHeight}	${window.innerWidth/window.innerHeight}`);	// Window is affected by dev. tools window
	canvas.width = window.innerWidth
	canvas.height = window.innerHeight
	
	document.body.appendChild(canvas)
	gl = canvas.getContext('webgl2')
	gl.clearColor(1, 0, 1, 1)
	
	console.log(gl.getParameter(gl.SHADING_LANGUAGE_VERSION));
	console.log(gl.getParameter(gl.VERSION));
	
	let vertexShader = gl.createShader(gl.VERTEX_SHADER)
	gl.shaderSource(vertexShader, vertShaderSource)
	gl.compileShader(vertexShader)

	if (!gl.getShaderParameter(vertexShader, gl.COMPILE_STATUS)) {
		console.log("vert: ", gl.getShaderInfoLog(vertexShader))
	}

	let fragmentShader = gl.createShader(gl.FRAGMENT_SHADER)
	gl.shaderSource(fragmentShader,fragShaderSource)
	gl.compileShader(fragmentShader)
	//console.log("frag.glsl: ", fragShaderSource)
	
	

	if (!gl.getShaderParameter(fragmentShader, gl.COMPILE_STATUS)) {
		console.log("frag: ", gl.getShaderInfoLog(fragmentShader))
	}

	program = gl.createProgram()
	gl.attachShader(program, vertexShader)
	gl.attachShader(program, fragmentShader)
	gl.linkProgram(program)

	let buffer = gl.createBuffer()
	gl.bindBuffer(gl.ARRAY_BUFFER, buffer)
	gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW)

	gl.useProgram(program)

	SetShaderUniforms()

	// Set shader attributes
	program.position = gl.getAttribLocation(program, 'pos')
	gl.enableVertexAttribArray(program.position)
	gl.vertexAttribPointer(program.position, 2, gl.FLOAT, false, 0, 0)

	draw()
}

function SetShaderUniforms() {
	program.resolution = gl.getUniformLocation(program, 'resolution')
	gl.uniform2fv(program.resolution, [canvas.width, canvas.height])
}

function draw() {
	// Reset Window to clear color(magenta)
	gl.clear(gl.COLOR_BUFFER_BIT)
	// Draw call
	gl.drawArrays(gl.TRIANGLES, 0, vertices.length / 2)

	
	//let aspectRatio = canvas.width/canvas.height
	//ctx.fillStyle = `rgb(${255},${255},${255})`
	//ctx.fillRect(points[0]/(aspectRatio > 1 ? 1 : aspectRatio)*canvas.width, points[1]*(aspectRatio > 1 ? aspectRatio : 1)*canvas.height, 20, 20)
}