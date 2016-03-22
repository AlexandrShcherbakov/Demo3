#include "SimpleInCube.h"

#include <iostream>

template <typename T>
uint argMax(T begin, T end) {
    uint maxIndex = 0;
    auto maxValue = begin;
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

std::tuple<std::vector<vec4>, std::vector<std::set<uint> >, std::vector<std::set<uint> > > createAdjacencyList(std::vector<Polygon>& poly) {
    std::vector<vec4> uniquePoints;
    std::vector<std::set<uint> > adjList;
    std::vector<std::set<uint> > owners;
    for (uint p = 0; p < poly.size(); ++p) {
        auto edges = poly[p].getEdges();
        for (uint edIt = 0; edIt < edges.size(); edIt += 2) {
            bool pointExist = false;
            uint indexA;
            for (indexA = 0; indexA < uniquePoints.size() && !pointExist; ++indexA)
                if (pointExist = (uniquePoints[indexA] == edges[edIt])) break;
            if (!pointExist) {
                uniquePoints.push_back(edges[edIt]);
                adjList.resize(uniquePoints.size());
                owners.resize(uniquePoints.size());
            }
            pointExist = false;
            uint indexB;
            for (indexB = 0; indexB < uniquePoints.size() && !pointExist; ++indexB)
                if (pointExist = (uniquePoints[indexB] == edges[edIt + 1])) break;
            if (!pointExist) {
                uniquePoints.push_back(edges[edIt + 1]);
                adjList.resize(uniquePoints.size());
                owners.resize(uniquePoints.size());
            }
            adjList[indexA].insert(indexB);
            adjList[indexB].insert(indexA);
            owners[indexA].insert(p);
            owners[indexB].insert(p);
        }
    }
    return std::make_tuple(uniquePoints, adjList, owners);
}

bool isBoundary(vec4 p, vec3 top, vec3 bottom) {
    return std::abs(p.x - top.x) < VEC_EPS ||
		   std::abs(p.y - top.y) < VEC_EPS ||
		   std::abs(p.z - top.z) < VEC_EPS ||
		   std::abs(p.x - bottom.x) < VEC_EPS ||
		   std::abs(p.y - bottom.y) < VEC_EPS ||
		   std::abs(p.z - bottom.z) < VEC_EPS;
}

std::vector<Polygon> simplifySmallCube(std::vector<Polygon>& poly,
                                     vec3 bottom,
                                     vec3 top) {
    std::vector<Polygon> newPoly(poly.begin(), poly.end());
    uint notBoundaryPoints;
    ///Create adjacency list
    auto adjInfo = createAdjacencyList(poly);
    auto points = std::get<0>(adjInfo);
    auto adjList = std::get<1>(adjInfo);
    auto owners = std::get<2>(adjInfo);
    ///Marker boundary points
    std::vector<uint> markers(points.size(), 0);///0 - unmarkered point
    notBoundaryPoints = points.size();
    ///1 - boundary point
    ///2 - deleted point
    std::queue<uint> queueToRemove;
    std::vector<bool> used(points.size(), false);
    std::vector<bool> workedPoly(poly.size(), false);
    for (uint i = 0; i < points.size(); ++i) {
        notBoundaryPoints -= (markers[i] = isBoundary(points[i], top, bottom));
        used[i] = markers[i];
        if (markers[i])
			for (auto p: owners[i]) {
                workedPoly[p] = true;
			}
    }
    for (uint i = 0; i < points.size(); ++i) {
        if (used[i])
            for (auto iter = adjList[i].begin(); iter != adjList[i].end(); ++iter)
                if (!used[*iter]) {
                    queueToRemove.push(*iter);
                    used[*iter] = true;
                }
    }
    ///While unmarkered points exists, marker its as deleted and merge with nearest boundary point
    //std::cout << notBoundaryPoints << std::endl;
    uint iterCount = notBoundaryPoints;
    while (queueToRemove.size() > 0) {
        /*if (100 * (iterCount - notBoundaryPoints + 1) / iterCount > 100 * (iterCount - notBoundaryPoints) / iterCount)
            std::cout << 100 * (iterCount - notBoundaryPoints) / iterCount << "% points removed" << std::endl;
        if (100 * (iterCount - notBoundaryPoints + 1) / iterCount > 100 * (iterCount - notBoundaryPoints) / iterCount)
            std::cout << queueToRemove.size() << std::endl;*/
        if (queueToRemove.size() == 0)
            break;
        uint target = queueToRemove.front();
        queueToRemove.pop();
        uint boundary = target;
        for (auto iter = adjList[target].begin(); iter != adjList[target].end(); ++iter) {
            if (markers[*iter] == 1) {
				if (boundary == target || length(points[boundary] - points[target]) > length(points[*iter] - points[target]))
					boundary = *iter;
            }
        }
		for (auto iter = adjList[target].begin(); iter != adjList[target].end(); ++iter)
            if (*iter != boundary)
                adjList[*iter].insert(boundary);
        for (auto p: owners[target]) {
            newPoly[p].changePoint(points[target], points[boundary]);
            workedPoly[p] = true;
        }
        markers[target] = 2;
        for (auto iter = adjList[target].begin(); iter != adjList[target].end(); ++iter)
            if (!used[*iter]) {
                queueToRemove.push(*iter);
                used[*iter] = true;
            }
		notBoundaryPoints--;
    }
    //std::cout << "OK" << std::endl;
    ///Create result
    std::vector<Polygon> result;
    for (uint p = 0; p < newPoly.size(); ++p)
        if (workedPoly[p] && !newPoly[p].isEmpty())
            result.push_back(newPoly[p]);
    //std::cout << "OK" << std::endl;
    return result;
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

std::vector<vec4> SimpleInCube(std::vector<vec4>& triangles, const uint maxDeep) {
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
