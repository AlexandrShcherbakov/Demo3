#include "SimpleInCube.h"

template <typename Iter>
uint argMax(Iter begin, Iter end) {
    uint maxIndex = 0;
    Iter maxValue = begin;
    uint count = 0;
    while (begin != end) {
        if (*maxValue < *begin) {
            maxValue = begin;
            maxIndex = count;
        }
        begin++;
        count++;
    }
    return maxIndex;
}

std::vector<Polygon> simplifyInSubcube(std::vector<Polygon>& poly,
                                     vec3 bottom,
                                     vec3 top,
                                     uint deep,
                                     uint maxDeep);

std::vector<Polygon> splitCube(std::vector<Polygon>& poly,
                                     vec3 bottom,
                                     vec3 top,
                                     uint deep,
                                     uint maxDeep) {
    float sides[] = {top.x - bottom.x, top.y - bottom.y, top.z - bottom.z};
    vec3 midTop = top;
    vec3 midBot = bottom;
    uint maxSide = argMax(sides, sides + 3);
    vec4 plane;
    if (maxSide == 0) {
        plane = vec4(1.0f, 0.0f, 0.0f, -(top.x + bottom.x) / 2.0f);
        midBot.x = midTop.x = (top.x + bottom.x) / 2.0f;
    } else if (maxSide == 1) {
        plane = vec4(0.0f, 1.0f, 0.0f, -(top.y + bottom.y) / 2.0f);
        midBot.y = midTop.y = (top.y + bottom.y) / 2.0f;
    } else {
        plane = vec4(0.0f, 0.0f, 1.0f, -(top.z + bottom.z) / 2.0f);
        midBot.z = midTop.z = (top.z + bottom.z) / 2.0f;
    }
    std::vector<Polygon> neg, pos;
    for (uint i = 0; i < poly.size(); ++i) {
        auto splited = poly[i].splitByPlane(plane);
        if (!splited.first.isEmpty())
            neg.push_back(splited.first);
        if (!splited.second.isEmpty())
            pos.push_back(splited.second);
    }
    std::vector<Polygon> result;
    if (pos.size() > 0)
        for (auto i: simplifyInSubcube(pos, midBot, top, deep + 1, maxDeep))
            result.push_back(i);
    if (neg.size() > 0)
        for (auto i: simplifyInSubcube(neg, bottom, midTop, deep + 1, maxDeep))
            result.push_back(i);
    return result;
}

std::vector<Polygon> simplifySmallCube(std::vector<Polygon>& poly,
                                     vec3 bottom,
                                     vec3 top) {
    ///TODO: fill this function
    return poly;
}

std::vector<Polygon> simplifyInSubcube(std::vector<Polygon>& poly,
                                     vec3 bottom,
                                     vec3 top,
                                     uint deep,
                                     uint maxDeep) {
    if (poly.size() == 0) return std::vector<Polygon>();
    if (deep < maxDeep) {
        return splitCube(poly, bottom, top, deep, maxDeep);
    } else {
        return simplifySmallCube(poly, bottom, top);
    }
}

std::vector<vec4> SimpleInCube(const std::vector<vec4>& triangles, const uint maxDeep) {
    std::vector<Polygon> poly;
    vec3 bottom(1.0f / VEC_EPS, 1.0f / VEC_EPS, 1.0f / VEC_EPS);
    vec3 top(-1.0f / VEC_EPS, -1.0f / VEC_EPS, -1.0f / VEC_EPS);
    for (uint i = 0; i < triangles.size(); i += 3) {
        for (uint j = 0; j < 3; ++j) {
            bottom.x = std::min(bottom.x, triangles[i + j].x);
            bottom.y = std::min(bottom.y, triangles[i + j].y);
            bottom.z = std::min(bottom.z, triangles[i + j].z);
            top.x = std::max(top.x, triangles[i + j].x);
            top.y = std::max(top.y, triangles[i + j].y);
            top.z = std::max(top.z, triangles[i + j].z);
        }
        poly.push_back(Polygon(triangles.begin() + i, triangles.begin() + i + 3));
    }
    std::vector<vec4> result;
    for (auto p: simplifyInSubcube(poly, bottom, top, 1, maxDeep))
        for (auto v: p.getEdges())
            result.push_back(v);
    return result;
}
