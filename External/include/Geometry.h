#ifndef GEOMETRY_H_INCLUDED
#define GEOMETRY_H_INCLUDED

#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <algorithm>
#include <fstream>

#include "VectorMath.h"

///*****************************************************************************
///Geometry class
///*****************************************************************************
/*
Набросок класса

Вектор4 Точки (множество)

Рёбра (множество из 2 вершин)
Треугольники (множество из 3 вершин)

Убрать внутренние вершины
Убрать маленькие треугольники
Конструктор из вектора вершин и вектора индексов

Слить 2 вершины

*/

class Geometry {
private:
    std::vector<vec4> points;
    //std::vector<uint> indices;
    //std::vector<uint> edges;
    std::vector<uint> triangles;
    //std::vector<uint> removedPoints;

    //void compressData();
    void compressIndices();
    void removeNullEdges();
    //void mergeSimilarPoints();
    void removeBadTriangles();
    void removeEqualEdges();
    vec3 positiveOrient(vec3 v);
    std::vector<std::vector<std::vector<uint> > > groupTrianglesByPlanes();

public:

    ///Small changes
    void mergeTwoPoints(uint p1, uint p2, vec4 newP); //OK
    void mergeTwoPoints(uint p1, uint p2); //OK
    void removeTriangle(uint p1, uint p2, uint p3); //OK

    ///Big changes
    void mergeSimilarPoints(); //OK
    void removePointsIntoPolygons();
    void removePointsIntoLines();
    void removeSmallTriangles(const float minSquare);
    void compressData();
    void removeExcessPoints();
    void removeInternalPoints();
    void FirstPlane();

    ///Loaders
    void loadFromTriangles(const std::vector<vec4> &points, const std::vector<uint>& indices); //OK
    void loadFromEdges(const std::vector<vec4> &edges);

	///Representations of geometry
	std::vector<vec4> getEdgeList();
	std::vector<uint> getEdgeListIndices();
	std::vector<std::set<uint> > getAdjacencyListIndices();
	std::vector<std::set<uint> > getBidirAdjacencyListIndices();

    ///Subsidary structure getters
    std::vector<std::vector<vec4> > getOutLines();

    ///Getters
    std::vector<vec4>& getPoints();
    std::vector<uint>& getTriangles();
    std::vector<uint>& getEdges();
    const std::vector<vec4>& getPoints() const;
    const std::vector<uint>& getTriangles() const;
    const std::vector<uint>& getEdges() const;
    std::vector<vec4> getEdgesAsPoints() const;
    void repairGeometry();

    uint getTrianglesNumber() const;
    uint getEdgesNumber() const;
    uint getPointsNumber() const;
    vec4 getTrianglePlant(const uint index) const;

    ///IO
    void Read(const std::string &path);
    void Write(const std::string &path) const;
};

#endif // GEOMETRY_H_INCLUDED
