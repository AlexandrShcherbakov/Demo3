#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include <utility>

#include "VectorMath.h"

class Polygon
{
    public:
        Polygon();
        template <typename Iter>
        Polygon(Iter begin, Iter end);
        Polygon(const Polygon &p);

        ///Tools
        std::pair<Polygon, Polygon> splitByPlane(const vec4& plane) const;
        bool isEmpty() const;

        ///Getters
        uint getSize() const;
        float getSquare() const;
        std::vector<vec4> getPoints() const;
        std::vector<vec4> getEdges() const;

    protected:
    private:
        std::vector<vec4> points;
};

#endif // POLYGON_H
