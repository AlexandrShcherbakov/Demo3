#include "Polygon.h"

Polygon::Polygon() {

}

template <class T>
Polygon::Polygon(T begin, T end) {
    points = std::vector<vec4>(begin, end);
}

Polygon::Polygon(const Polygon &p) {
    points.clear();
    for (auto i: p.points) {
        points.push_back(i);
    }
}

uint Polygon::getSize() const {
    return points.size();
}

float Polygon::getSquare() const {
    float square = 0.0f;
    for (uint i = 2; i < points.size(); ++i) {
        square += length(cross((points[i] - points[0]).xyz(), (points[i - 1] - points[0]).xyz()));
    }
    return std::abs(square) / 2.0f;
}

std::vector<vec4> Polygon::getPoints() const {
    return points;
}

bool Polygon::isEmpty() const {
    return points.size() < 3;
}

std::pair<Polygon, Polygon> Polygon::splitByPlane(const vec4& plane) const {
    ///If there are no points in polygon
    if (getSize() == 0) return std::make_pair(Polygon(), Polygon());

    std::vector<vec4> newPoints;

    for (uint i = 1; i <= points.size(); ++i) {
        newPoints.push_back(points[i - 1]);
        float t = -dot(points[i - 1], plane) / dot((points[i % points.size()] - points[i - 1]).xyz(), plane.xyz());
        if (t > 0 && t < 1) {
            newPoints.push_back(points[i - 1] + (points[i % points.size()] - points[i - 1]) * t);
            newPoints.back().w = 1.0f;
        }
    }

    std::vector<vec4> neg, pos;
    for (auto p: newPoints) {
        if (dot(p, plane) <= VEC_EPS)
            neg.push_back(p);
        if (dot(p, plane) >= -VEC_EPS)
            pos.push_back(p);
    }
    return std::make_pair(Polygon(neg.begin(), neg.end()), Polygon(pos.begin(), pos.end()));
}

std::vector<vec4> Polygon::getEdges() const {
    std::vector<vec4> edges;
    if (points.size() < 2) return edges;
    edges.push_back(points.back());
    for (uint i = 0; i < points.size() - 1; ++i) {
        edges.push_back(points[i]);
        edges.push_back(points[i]);
    }
    if (edges.size() < 3) return edges;
    edges.push_back(points.back());
    return edges;
}

void Polygon::changePoint(vec4 oldValue, vec4 newValue) {
    for(uint i = 0; i < points.size(); ++i) {
		if (points[i] == oldValue) {
            points[i] = newValue;
            if (points[(points.size() + i - 1) % points.size()] == points[i]) i = (points.size() + i - 1) % points.size();
            if (points[(i + 1) % points.size()] == points[i]) {
                points.erase(points.begin() + i);
                i--;
            }
		}
    }
}
