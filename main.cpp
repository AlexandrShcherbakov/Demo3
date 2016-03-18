#include "HydraExport.h"
#include <iostream>
#include <fstream>
#include <string>
#include "VectorMath.h"
#include "Simplifier.h"
#include "Geometry.h"
#include "SceneGeometry.h"
#include "GL\glew.h"
#include "GLUtility.h"
#include "GL\freeglut.h"

using namespace std;

HydraGeomData data;
Geometry baseGeom;
Geometry simple;
SceneGeometry g1, g2;
mat4 camera;
vec3 rot;
ShaderProgram *p1, *p2;

int flag = 1;

GLuint Prog1;
GLuint vao;
GLuint pnt, ind;

bool startSiplify = false;

void RenderLayouts() {
    glClear(GL_COLOR_BUFFER_BIT);
    p1->setUniform("cam", camera * Matrix4Rotate(rot));
    if (flag & 1)
    p1->Draw();
    if (startSiplify) {
		g2.removeSmallestEdges(100);
		p2->setGeometry(g2);
    }
    p2->setUniform("cam", camera * Matrix4Rotate(rot));
    if (flag & 2)
    p2->Draw();
	glutSwapBuffers();
}

void KeyboardEvents(unsigned char key, int x, int y) {
	if (key == 'w' || key == 'W') {
        camera *= Matrix4Shift(vec3(0.0f, 0.0f, 0.01f));
	} else if (key == 's' || key == 'S') {
		camera *= Matrix4Shift(vec3(0.0f, 0.0f, -0.01f));
	} else if (key == 'a' || key == 'A') {
        camera *= Matrix4Shift(vec3(0.01f, 0.0f, 0.0f));
	} else if (key == 'd' || key == 'D') {
	    camera *= Matrix4Shift(vec3(-0.01f, 0.0f, 0.0f));
	} else if (key == '1') {
        flag = 1;
	} else if (key == '2') {
        flag = 2;
	} else if (key == '3') {
        flag = 3;
	} else if (key == '+') {
		g2.removeSmallestEdges(100);
		p2->setGeometry(g2);
		startSiplify = true;
	}
}

void SpecialButtons(int key, int x, int y) {
    if (key == GLUT_KEY_RIGHT) {
        rot.y += 0.1f;
	} else if (key == GLUT_KEY_LEFT) {
	    rot.y -= 0.1f;
	} else if (key == GLUT_KEY_UP) {
	    rot.x += 0.1f;
	} else if (key == GLUT_KEY_DOWN) {
        rot.x -= 0.1f;
	}
}

void IdleFunc() {
    glutPostRedisplay();
}

void MouseMove(int x, int y) {
}

void MouseClick(int button, int state, int x, int y) {
}

void InitializeGLUT(int argc, char **argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitContextVersion(3, 0);
	glutInitWindowPosition(-1, -1);
	glutInitWindowSize(800, 600);
	glutCreateWindow("Demo 3");

	glutDisplayFunc(RenderLayouts);
    glutKeyboardFunc(KeyboardEvents);
    glutSpecialFunc(SpecialButtons);
    glutIdleFunc(IdleFunc);
    glutMotionFunc(MouseMove);
    glutMouseFunc(MouseClick);
}

SceneGeometry ExtractSceneFromFile(const std::string &path) {
	data.read(path);
	std::cout << data.getVerticesNumber() << ' ' << data.getIndicesNumber() << "\n";
 	std::vector<vec4> points;
	for (uint i = 0; i < data.getVerticesNumber(); i++) {
		const float *v = data.getVertexPositionsFloat4Array() + i * 4;
		points.push_back(vec4(v[0], v[1], v[2], 1));
	}

	std::vector<uint> indices(data.getIndicesNumber());
	for (uint i = 0; i < indices.size(); ++i)
		indices[i] = data.getTriangleVertexIndicesArray()[i];
    SceneGeometry baseGeom;
    //points.resize(600);
    //indices.resize(600);
    baseGeom.setTriangles(points, indices);
    simple.loadFromTriangles(points, indices);
    return baseGeom;
}

Geometry ReduceTriangles(const Geometry &geom) {
    Geometry res = geom;
    std::vector<vec4> rfTri = res.getPoints();
    for (uint i = 0; i < res.getTrianglesNumber(); ++i) {
        vec4 mid = (rfTri[3 * i] + rfTri[3 * i + 1] + rfTri[3 * i + 2]) / 3;
        rfTri[3 * i + 0] = (rfTri[3 * i + 0] + mid) / 2;
        rfTri[3 * i + 1] = (rfTri[3 * i + 1] + mid) / 2;
        rfTri[3 * i + 2] = (rfTri[3 * i + 2] + mid) / 2;
    }
    return res;
}

void InitializeGL() {
	glewInit();
    //Prog1 = LoadShader("GreenSkeleton");
    p1 = new ShaderProgram("GreenSkeleton");
    p2 = new ShaderProgram("GreenSkeleton");
    Matrix4Perspective(camera, 90, 1, 0.05, 50);
    camera *= Matrix4Shift(vec3(0.0f, 0.0f, -1.0f));
}

void CreateBuffersFromData() {
	pnt = CreateGLBuffer("points", GL_ARRAY_BUFFER, data.getVerticesNumber() * 4 * sizeof(float), data.getVertexPositionsFloat4Array());
	ind = CreateGLBuffer("indices", GL_ELEMENT_ARRAY_BUFFER, data.getIndicesNumber() * sizeof(uint), data.getTriangleVertexIndicesArray());
	glGenVertexArrays(1, &vao);                                 CHECK_GL_ERRORS
    glBindVertexArray(vao);                                     CHECK_GL_ERRORS
    auto pointsLocation = glGetAttribLocation(3, "points"); CHECK_GL_ERRORS
    glBindBuffer(GL_ARRAY_BUFFER, pnt);                                CHECK_GL_ERRORS
    glEnableVertexAttribArray(pointsLocation);                                   CHECK_GL_ERRORS
    glVertexAttribPointer(pointsLocation, 4, GL_FLOAT, GL_FALSE, 0, 0);          CHECK_GL_ERRORS

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ind);                       CHECK_GL_ERRORS
    glBindVertexArray(vao);
}


int main(int argc, char **argv) {
	InitializeGLUT(argc, argv);
	//std::cout << "OK" << std::endl;
	InitializeGL();
	//baseGeom = ExtractSceneFromFile("Scenes/dabrovic-sponza/scene.vsgf");
	g1 = ExtractSceneFromFile("Scenes/dabrovic-sponza/scene.vsgf");
	//baseGeom.Write("Step 1 - unique points");
	//baseGeom = ReduceTriangles(baseGeom);
	//CreateBuffersFromData();
	//baseGeom.repairGeometry();
	//baseGeom.Write("Step 1 - unique points");
    //baseGeom.Read("Step 1 - unique points");
    //return 0;
	//baseGeom = ReadGeometry("Step 1 - unique points");
	//baseGeom = NewQuadricErrorMetricsSimplifier(geom, 1);
	//WriteGeometry(baseGeom, "Step 4 - remove bad edges");
	//baseGeom = ReadGeometry("Step 1 - unique points");
	//baseGeom.repairGeometry();
	//simple = ReadGeometry("Step 4 - remove bad edges");
	//simple = NewQuadricErrorMetricsSimplifier(baseGeom, 3000);
	//simple = SmallEdgesSimplifier(baseGeom, 300);
	//simple.mergeSimilarPoints();
	//simple.Write("Step 1 - unique points");
	//simple.removeInternalPoints();
	//baseGeom.FirstPlane();
    g2 = g1;
    g2.removeSimilarPoints();
    std::cout << g2.getEdges().size() << " edges" << std::endl;
    //g2.removeSmallestEdges(90000);
    std::cout << g2.getEdges().size() << " edges" << std::endl;
	p1->setGeometry(g1);
    p1->setUniform("col", vec3(1, 0, 0));
    p2->setGeometry(g2);
    p2->setUniform("col", vec3(0, 1, 0));
	glutMainLoop();
    //Geometry gm = NewQuadricErrorMetricsSimplifier(baseGeom, 3000);
}
