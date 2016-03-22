#include "SceneGeometry.h"

#include <iostream>

std::unordered_map<uint, std::set<uint> > CreateAdjacencyList(const std::vector<uint>& triangles) {
    std::unordered_map<uint, std::set<uint> > adjacencyList;
    for (uint i = 0; i < triangles.size(); i += 3) {
        adjacencyList[triangles[i + 0]].insert(triangles[i + 1]);
        adjacencyList[triangles[i + 0]].insert(triangles[i + 2]);
        adjacencyList[triangles[i + 1]].insert(triangles[i + 2]);
    }
    return adjacencyList;
}


void SceneGeometry::setTriangles(const std::vector<vec4>& points, const std::vector<uint>& triangles) {
	///Clear internal containers
    this->points.clear();
    this->edges.clear();
    this->edgesOfPoints.clear();
    edgeKey = 0;

    ///Add points to container
    for (uint i = 0; i < points.size(); ++i)
        this->points[i] = points[i];

	///Add edges and relate it with points
    auto adjacencyList = CreateAdjacencyList(triangles);
    for (auto i = adjacencyList.begin(); i != adjacencyList.end(); ++i) {
        uint vert1 = i->first;
        for (auto j = adjacencyList[vert1].begin(); j != adjacencyList[vert1].end(); ++j) {
            std::vector<uint> edge = {vert1, *j};
            uint index = edgeKey = edges.size();
            edges[index] = edge;
            edgesOfPoints.emplace(vert1, index);
            edgesOfPoints.emplace(*j, index);
        }
    }
}

void SceneGeometry::setEdges(const std::vector<vec4>& points, const std::vector<uint>& edges) {
	///Clear internal containers
    this->points.clear();
    this->edges.clear();
    this->edgesOfPoints.clear();
    edgeKey = 0;

    ///Add points to container
    for (uint i = 0; i < points.size(); ++i)
        this->points[i] = points[i];

	///Add edges and relate it with points
    for (uint i = 0; i < edges.size(); i += 2) {
		std::vector<uint> edge = {edges[i], edges[i + 1]};
		uint index = edgeKey = this->edges.size();
		this->edges[index] = edge;
		edgesOfPoints.emplace(edges[i + 0], index);
		edgesOfPoints.emplace(edges[i + 1], index);
    }
}

std::vector<vec4> SceneGeometry::getPoints() const {
    std::vector<vec4> points;

    for (auto i = this->points.begin(); i != this->points.end(); ++i) {
		///Update size
        if (points.size() <= i->first) points.resize(i->first + 1);
		///Add point
        points[i->first] = i->second;
    }
    return points;
}

std::vector<uint> SceneGeometry::getEdges() const {
    std::vector<uint> edges;
    for (auto i = this->edges.begin(); i != this->edges.end(); ++i) {
        edges.push_back(i->second[0]);
        edges.push_back(i->second[1]);
    }
    return edges;
}

void SceneGeometry::removeSimilarPoints() {
	std::vector<vec4> uniquePoints;
	std::vector<uint> firstEntry;

	std::vector<vec4> points = getPoints();

    for (uint i = 0; i < points.size(); ++i) {
        if (100 * i / points.size() > 100 * (i - 1) / points.size())
			std::cout << 100 * i / points.size() << "% points updated\n";
		bool isUnique = true;
        for (uint j = uniquePoints.size(); j > 0 && isUnique; --j) {
            isUnique = this->points[i] != uniquePoints[j - 1];
            ///If this point already found
			if (!isUnique)
				mergeTwoPoints(i, firstEntry[j - 1]);
        }
        ///If it is first entry of this point
        if (isUnique) {
            firstEntry.push_back(i);
            uniquePoints.push_back(this->points[i]);
        }
    }
}

void SceneGeometry::mergeTwoPoints(uint p1, uint p2) {
    mergeTwoPoints(p1, p2, (points[p1] + points[p2]) / 2);
}

void SceneGeometry::mergeTwoPoints(uint p1, uint p2, vec4 v) {
    std::map<uint, std::pair<uint, uint> > adjacencyList;

    ///Sort argument indices
    if (p1 > p2) std::swap(p1, p2);

	///Choose new edges
    std::map<uint, std::vector<uint> > newEdges;
    auto iter = edgesOfPoints.find(p1);
    for (uint i = 0; i < edgesOfPoints.count(p1); ++i, ++iter) {
        uint edgeKey = iter->second;
    	if (edges[edgeKey][0] == p2) edges[edgeKey][0] = p1;
        if (edges[edgeKey][1] == p2) edges[edgeKey][1] = p1;
		if (edges[edgeKey][0] == edges[edgeKey][1]) continue;
        if (edges[edgeKey][0] > edges[edgeKey][1]) std::swap(edges[edgeKey][0], edges[edgeKey][1]);
        bool isNew = true;
        for (uint j = 0; j < 2; ++j)
			if (edges[edgeKey][j] != p1)  {
				auto tarIter = edgesOfPoints.find(edges[edgeKey][j]);
				while (tarIter->second != edgeKey) tarIter++;
				edgesOfPoints.erase(tarIter);
			}
        for (auto j = newEdges.begin(); j != newEdges.end() && isNew; ++j)
            isNew = edges[edgeKey][0] != j->second[0] || edges[edgeKey][1] != j->second[1];
		if (isNew)
            newEdges[edgeKey] = edges[edgeKey];
    }
    iter = edgesOfPoints.find(p2);
    for (uint i = 0; i < edgesOfPoints.count(p2); ++i, ++iter) {
    	uint edgeKey = iter->second;
    	if (edges[edgeKey][0] == p2) edges[edgeKey][0] = p1;
        if (edges[edgeKey][1] == p2) edges[edgeKey][1] = p1;
		if (edges[edgeKey][0] == edges[edgeKey][1]) continue;
        if (edges[edgeKey][0] > edges[edgeKey][1]) std::swap(edges[edgeKey][0], edges[edgeKey][1]);
		for (uint j = 0; j < 2; ++j)
			if (edges[edgeKey][j] != p1)  {
				auto tarIter = edgesOfPoints.find(edges[edgeKey][j]);
				while (tarIter->second != edgeKey) tarIter++;
				edgesOfPoints.erase(tarIter);
			}
        bool isNew = true;
        for (auto j = newEdges.begin(); j != newEdges.end() && isNew; ++j)
            isNew = edges[edgeKey][0] != j->second[0] || edges[edgeKey][1] != j->second[1];
		if (isNew)
            newEdges[edgeKey] = edges[edgeKey];
    }

	///Update edges
	iter = edgesOfPoints.find(p1);
    for (uint i = 0; i < edgesOfPoints.count(p1); ++i, ++iter)
        edges.erase(iter->second);
    iter = edgesOfPoints.find(p2);
    for (uint i = 0; i < edgesOfPoints.count(p2); ++i, ++iter)
    	edges.erase(iter->second);
	edgesOfPoints.erase(p2);
	edgesOfPoints.erase(p1);
    for (auto i = newEdges.begin(); i != newEdges.end(); ++i) {
        edges[i->first] = i->second;
        edgesOfPoints.emplace(i->second[0], i->first);
        edgesOfPoints.emplace(i->second[1], i->first);
    }

	///Remove second point
    points.erase(p2);

	///Update first point
    points[p1] = v;
}

void SceneGeometry::removeSmallestEdgesOneVector(const uint count) {
	///Create vector of target edges
    std::vector<std::pair<float, uint> > edges;
	for (auto i = this->edges.begin(); i != this->edges.end(); ++i) {
		float l1 = length(points[i->second[0]] - points[i->second[1]]);
		edges.push_back(std::make_pair(l1, i->first));
	}
    sort(edges.begin(), edges.end());
    edges.resize(5 * count);

    uint was = this->edges.size();

    ///Iteratively remove edges
    for (uint iter = 0; was - this->edges.size() < count; iter++) {
		///Merge points
        mergeTwoPoints(this->edges[edges[0].second][0], this->edges[edges[0].second][1]);

		///Update vector
		std::vector<std::pair<float, uint> > newEdges;
        for (uint i = 1; i < edges.size(); ++i) {
            if (this->edges.find(edges[i].second) == this->edges.end()) continue;
            float l1 = length(points[this->edges[edges[i].second][0]] - points[this->edges[edges[i].second][1]]);
            uint index = newEdges.size();
            newEdges.push_back(std::make_pair(l1, edges[i].second));
            while (index > 0 && newEdges[index - 1].first > newEdges[index].first)
				std::swap(newEdges[index - 1], newEdges[index]);
        }
        edges = newEdges;
    }
}


void SceneGeometry::removeSmallestEdges(uint count = 1) {
    if (count < 1000) removeSmallestEdgesOneVector(count);
    else {
		double v = edges.size() * (1 + log2(edges.size()));
        uint k = (uint)round(pow(10 * count * count / v, 1 / 3.0));
		std::cout << k << std::endl;
        for (uint i = 0; i < k; ++i) {
			if (100 * i / k > 100 * (i - 1) / k)
				std::cout << 100 * i / k << "% of edges removed " << count / k << std::endl;
            removeSmallestEdgesOneVector(count / k);
			std::cout << edges.size() << std::endl;
        }
    }

	/*///sort edges
    std::vector<std::pair<float, std::pair<uint, uint> > > edges;
	for (auto i = this->edges.begin(); i != this->edges.end(); ++i) {
		float l1 = length(points[i->second[0]] - points[i->second[1]]);
		edges.push_back(std::make_pair(l1, std::make_pair(i->second[0], i->second[1])));
	}
    sort(edges.begin(), edges.end());

	///Remove edges
    std::set<uint> removedPoints;
    for (uint i = 0; i < count; ++i) {
        if (removedPoints.find(edges[i].second.first) != removedPoints.end()) continue;
        if (removedPoints.find(edges[i].second.second) != removedPoints.end()) continue;
        removedPoints.insert(edges[i].second.second);
        mergeTwoPoints(edges[i].second.first, edges[i].second.second);
    }*/

    /*for (uint iter = 0; iter < count; ++iter) {
        if (100 * iter / count > 100 * (iter - 1) / count)
			std::cout << 100 * iter / count << "% edges removed" << std::endl;
        ///Find smallest edge
        float minLength = 1e9f;
        uint edgeIndex = 0;
        for (auto i = edges.begin(); i != edges.end(); ++i) {
            float l1 = length(points[i->second[0]] - points[i->second[1]]);
            if (minLength > l1) {
                minLength = l1;
                edgeIndex = i->first;
            }
        }

        ///Remove edge
        mergeTwoPoints(edges[edgeIndex][0], edges[edgeIndex][1]);
    }*/
}

void SceneGeometry::squeezeData() {
    std::vector<vec4> edges;
    for (auto i = this->edges.begin(); i != this->edges.end(); ++i) {
        edges.push_back(points[i->second[0]]);
        edges.push_back(points[i->second[1]]);
    }
    setEdges(edges);
}

void SceneGeometry::setEdges(const std::vector<vec4>& edges) {
	///Clear internal containers
    this->points.clear();
    this->edges.clear();
    this->edgesOfPoints.clear();
    edgeKey = 0;

    ///Set edges and create vector of points
    std::vector<vec4> uniquePoints;
    std::vector<uint> newEdges;
    for (uint i = 0; i < edges.size(); ++i) {
		bool firstEntry = true;
        uint index;
        for (uint j = uniquePoints.size(); j > 0 && firstEntry; --j) {
            if (!(firstEntry = (edges[i] != uniquePoints[j])))
                index = j;
        }
        if (firstEntry) {
            index = uniquePoints.size();
            uniquePoints.push_back(edges[i]);
        }
        newEdges.push_back(index);
    }
    std::cout << uniquePoints.size() << ' ' << newEdges.size() << std::endl;
    setEdges(uniquePoints, newEdges);
}
