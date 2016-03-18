#ifndef SCENEGEOMETRY_H
#define SCENEGEOMETRY_H

#include <vector>
#include <unordered_map>
#include <set>
#include <map>
#include <algorithm>

#include "VectorMath.h"

class SceneGeometry
{
	public:
		///Setters
		void setTriangles(const std::vector<vec4>& points, const std::vector<uint>& triangles);
		void setEdges(const std::vector<vec4>& edges);
		void setEdges(const std::vector<vec4>& points, const std::vector<uint>& edges);

		///Getters
        std::vector<vec4> getPoints() const;
        std::vector<uint> getEdges() const;

        ///Transformers
        void removeSimilarPoints();
        void removeSmallestEdges(uint count);
        void squeezeData();
	protected:
	private:
		uint edgeKey;

		///Containers
        std::unordered_map<uint, vec4> points;
        std::unordered_multimap<uint, uint> edgesOfPoints;
        std::unordered_map<uint, std::vector<uint> > edges;

		///Internal tools
        void mergeTwoPoints(uint p1, uint p2);
        void mergeTwoPoints(uint p1, uint p2, vec4 v);
        void removeSmallestEdgesOneVector(const uint count);
};

#endif // SCENEGEOMETRY_H
