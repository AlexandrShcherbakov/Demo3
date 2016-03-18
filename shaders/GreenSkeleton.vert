#version 330

in vec4 points;

uniform mat4 cam;

void main() {
	gl_Position = cam * points;
}
