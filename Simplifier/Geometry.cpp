#include "Geometry.h"

///Small changes
void Geometry::mergeTwoPoints(uint p1, uint p2) {
    mergeTwoPoints(p1, p2, (points[p1] + points[p2]) / 2.0f);
}

void Geometry::mergeTwoPoints(uint p1, uint p2, vec4 newP) {
    ///Update points
    /*points[p1] = newP;
    indices[p2] = p1;*/
	uint tmp = std::min(p1, p2);
	p2 = std::max(p1, p2);
	p1 = tmp;

	std::vector<uint> newTriangles;
    for (uint i = 0; i < triangles.size(); i += 3) {
        if (triangles[i + 0] == p2) triangles[i + 0] = p1;
        if (triangles[i + 1] == p2) triangles[i + 1] = p1;
        if (triangles[i + 2] == p2) triangles[i + 2] = p1;
        sort(triangles.begin() + i, triangles.begin() + i + 3);
        if (triangles[i + 0] == triangles[i + 1] || triangles[i + 1] == triangles[i + 2]) continue;
        newTriangles.push_back(triangles[i + 0]);
        newTriangles.push_back(triangles[i + 1]);
        newTriangles.push_back(triangles[i + 2]);
    }
    points.erase(points.begin() + p2);
    points[p1] = newP;
}


void Geometry::removeTriangle(uint p1, uint p2, uint p3) {
    vec4 resP = (points[p1] + points[p2] + points[p3]) / 3.0f;
    mergeTwoPoints(p1, p2, resP);
    mergeTwoPoints(p1, p3, resP);
}

void Geometry::removePointsIntoPolygons() {
    ///Create sets of plants for points and sets of edges
    std::vector<std::vector<vec4> > plants(points.size());
    std::vector<std::set<uint> > edgeList(points.size());
    for (uint i = 0; i < triangles.size(); i += 3) {
		sort(triangles.begin() + i, triangles.begin() + i + 3);
        uint a_ind = triangles[i + 0];
        uint b_ind = triangles[i + 1];
        uint c_ind = triangles[i + 2];
        vec4 plant = getTrianglePlant(i / 3);
        bool flag = false;
        for (uint i = 0; i < plants[a_ind].size() && !flag; ++i) {
            flag = flag || plant == plants[a_ind][i];
        }
        if (!flag) {
            plants[a_ind].push_back(plant);
        }
        flag = false;
        for (uint i = 0; i < plants[b_ind].size() && !flag; ++i) {
            flag = flag || plant == plants[b_ind][i];
        }
        if (!flag) {
            plants[b_ind].push_back(plant);
        }
        flag = false;
        for (uint i = 0; i < plants[c_ind].size() && !flag; ++i) {
            flag = flag || plant == plants[c_ind][i];
        }
        if (!flag) {
            plants[c_ind].push_back(plant);
        }
        edgeList[a_ind].insert(b_ind);
        edgeList[a_ind].insert(c_ind);
        edgeList[b_ind].insert(a_ind);
        edgeList[b_ind].insert(c_ind);
        edgeList[c_ind].insert(a_ind);
        edgeList[c_ind].insert(b_ind);
    }

    ///Remove bad points
    for (uint i = 0; i < plants.size(); ++i) {
		//std::cout << i << ' ' << plants[i].size() << ' ' << edgeList[i].size() << std::endl;
        if (100 * (i + 1) / plants.size()  > 100 * i / plants.size())
			std::cout << 100 * i / plants.size() << "% bad points removed!" << std::endl;
        if (plants[i].size() != 1) continue;
        auto j = edgeList[i].begin();
        uint mn = *j;
        for (j++; j != edgeList[i].end(); ++j) {
            if (length(points[i] - points[mn]) > length(points[i] - points[*j])) {
                mn = *j;
            }
        }
		//std::cout << "Nearest point found" << std::endl;
        mergeTwoPoints(mn, i, points[mn]);
        for (uint j = i + 1; j < edgeList.size(); ++j) {
            auto badNumIt = edgeList[j].find(i);
            if (badNumIt != edgeList[j].end()) {
				edgeList[j].erase(badNumIt);
                if (j != mn) edgeList[j].insert(mn);
				edgeList[j].insert(mn);
            }
        }
        //std::cout << i << " removed" << std::endl;
    }
    compressData();
    repairGeometry();
}

vec3 Geometry::positiveOrient(vec3 v) {
    if (v.x < 0) return v * -1;
    if (v.x > 0) return v;
    if (v.y < 0) return vec3(0, -v.y, -v.z);
    if (v.y > 0) return v;
    if (v.z < 0) return vec3(0, 0, -v.z);
    if (v.z > 0) return v;
    return v;
}

void Geometry::removePointsIntoLines() {
	std::cout << "Remove points into lines start: " << points.size() << ' ' << triangles.size() << std::endl;
    auto pointLines = getOutLines();
    auto edgeList = getBidirAdjacencyListIndices();
    std::vector<uint> IndexMap(points.size());
    for (uint i = 0; i < pointLines.size(); ++i) {
		IndexMap[i] = i;
        if (pointLines[i].size() == 1) {
        	int nearest = -1;
            for (auto j = edgeList[i].begin(); j != edgeList[i].end(); ++j) {
                if (IndexMap[*j] == i) continue;
                if (nearest == -1 || length(points[i] - points[*j]) < length(points[i] - points[nearest])) {
                    nearest = *j;
                }
            }
            IndexMap[i] = nearest;
        }
    }
    for (uint i = 0; i < IndexMap.size(); ++i) {
		uint j = i;
        while (j != IndexMap[j]) {
			j = IndexMap[j];
			std::cout << j << ' ' << IndexMap[j] << std::endl;
        }
		IndexMap[i] = j;
		points[i] = points[j];
    }
    mergeSimilarPoints();
    std::cout << "Remove points into lines end: " << points.size() << ' ' << triangles.size() << std::endl;
    /*///Create sets of lines for points and sets of edges
    //std::vector<std::set<vec3> > lines(points.size());
    std::vector<std::vector<vec3> > lines(points.size());
    std::vector<std::set<uint> > edgeList(points.size());
    for (uint i = 0; i < edges.size(); i += 2) {
        uint a_ind = edges[i + 0];
        uint b_ind = edges[i + 1];
        vec3 a = points[a_ind].xyz();
        vec3 b = points[b_ind].xyz();
        vec3 line = normalize(b - a);
        bool flag = false;
        for (uint j = 0; j < lines[a_ind].size() && !flag; ++j) {
            flag = flag || length(cross(lines[a_ind][j], line)) < VEC_EPS;
        }
        if (!flag) {
            lines[a_ind].push_back(line);
            edgeList[a_ind].insert(b_ind);
        }
        flag = false;
        for (uint j = 0; j < lines[b_ind].size() && !flag; ++j) {
            flag = flag || length(cross(lines[b_ind][j], line)) < VEC_EPS;
        }
        if (!flag) {
            lines[b_ind].push_back(line);
            edgeList[b_ind].insert(a_ind);
        }
    }

    ///Remove bad points
    for (uint i = 0; i < points.size(); ++i) {
        if (lines[i].size() != 1) continue;
        auto j = edgeList[i].begin();
        int mn = -1;
        for (; j != edgeList[i].end(); ++j) {
            //vec3 line = normalize(points[i].xyz() - points[*j].xyz());
            vec3 line = lines[i][0];
            if (mn == -1 || (length(points[i] - points[mn]) > length(points[i] - points[*j]) && length(cross(points[i].xyz() - points[mn].xyz(), line)) < VEC_EPS * 1e5)) {
                mn = *j;
            }
        }
        mergeTwoPoints(mn, i, points[mn]);
    }
    compressData();
    repairGeometry();*/
}

void Geometry::compressData() {
    ///Create new indices
    std::vector<vec4> newPoints;
    std::vector<uint> newIndices;
    std::vector<uint> indexMap(points.size());
    /*for (uint i = 0; i < points.size(); ++i) {
        if (indices[i] == i) {
            indexMap[i] = newIndices.size();
            newIndices.push_back(newPoints.size());
            newPoints.push_back(points[i]);
        }
    }
    for (uint i = 0; i < indices.size(); ++i) {
        if (indices[i] != i) {
            indexMap[i] = indexMap[indices[i]];
        }
    }
    for (uint i = 0; i < edges.size(); ++i) {
        edges[i] = indexMap[edges[i]];
    }
    for (uint i = 0; i < triangles.size(); ++i) {
        triangles[i] = indexMap[triangles[i]];
    }
    points = newPoints;
    indices = newIndices;
    ////Create index map
    //std::sort(removedPoints.begin(), removedPoints.end());
    std::vector<uint> newIndices(points.size());
    for (uint i = 0, sh = 0; i < points.size(); ++i) {
        if (removedPoints.size() > sh && i == removedPoints[sh]) sh++;
        newIndices[i] = i - sh;
    }

    ///Update edges
    std::vector<uint> newEdges;
    for (uint i = 0; i < edges.size(); i += 2) {
        if (edges[i] == edges[i + 1]) continue;
        newEdges.push_back(newIndices[edges[i]]);
        newEdges.push_back(newIndices[edges[i + 1]]);
    }
    edges = newEdges;

    ///Update triangles
    std::vector<uint> newTriangles;
    for (uint i = 0; i < triangles.size(); i += 3) {
        if (triangles[i] == triangles[i + 1] ||
            triangles[i] == triangles[i + 2] ||
            triangles[i + 1] == triangles[i + 2]) continue;
        newTriangles.push_back(newIndices[triangles[i]]);
        newTriangles.push_back(newIndices[triangles[i + 1]]);
        newTriangles.push_back(newIndices[triangles[i + 2]]);
    }
    triangles = newTriangles;

    ///Update points
    std::vector<vec4> newPoints;
    for (uint i = 0, sh = 0; i < points.size(); ++i) {
        if (removedPoints.size() > sh && i == removedPoints[sh]) sh++;
        else newPoints.push_back(points[i]);
    }
    points = newPoints;

    removedPoints.clear();*/
}

void Geometry::compressIndices() {
	/*std::vector<vec4> edgeBegin;
	std::vector<std::vector<vec4> > edgeEnds;
    for (uint i = 0; i < this->edges.size(); i += 2) {
        vec4 a;// = points[this->edges[i]];
		vec4 b;// = points[this->edges[i + 1]];
		if (a == b) continue;
        bool flag1 = false;
        for (uint j = 0; j < edgeBegin.size() && !flag1; ++j) {
            if (a == edgeBegin[j]) {
                flag1 = true;
                bool flag2 = false;
                for (uint h = 0; h < edgeEnds[j].size() && !flag2; ++h) {
                    flag2 = flag2 || edgeEnds[j][h] == b;
                }
                if (!flag2) edgeEnds[j].push_back(b);
            }
        }
    }
    std::vector<vec4> edges;
    for (auto i = 0; i < edgeBegin.size(); ++i) {
        for (auto j = 0; j < edgeEnds[i].size(); ++j) {
            edges.push_back(edgeBegin[i]);
            edges.push_back(edgeEnds[i][j]);
        }
    }
    loadFromEdges(edges);*/
}

void Geometry::loadFromTriangles(const std::vector<vec4> &points, const std::vector<uint>& indices) {
    this->points = std::vector<vec4>(points.begin(), points.end());
    this->triangles = std::vector<uint>(indices.begin(), indices.end());
    /*this->indices = std::vector<uint>(this->points.size());
    for (uint i = 0; i < this->indices.size(); ++i)
        this->indices[i] = i;
    edges.clear();
    std::vector<std::set<uint> > edgeList(points.size());
    for (uint i = 0; i < triangles.size(); i += 3) {
		sort(triangles.begin() + i, triangles.begin() + i + 3);
		uint a = triangles[i];
        uint b = triangles[i + 1];
        uint c = triangles[i + 2];
        if (a < b)
			edgeList[a].insert(b);
		if (a < c)
			edgeList[a].insert(c);
		if (b < c)
			edgeList[b].insert(c);
    }
    for (uint i = 0; i < edgeList.size(); ++i) {
        for (auto j = edgeList[i].begin(); j != edgeList[i].end(); ++j) {
            edges.push_back(i);
            edges.push_back(*j);
        }
    }*/
}

void Geometry::loadFromEdges(const std::vector<vec4> &edges) {
    std::vector<vec4> points;
    std::vector<uint> edgesIndices;
    for (uint i = 0; i < edges.size(); ++i) {
        bool flag = false;
        for (uint j = 0; j < points.size() && !flag; ++j) {
            flag = edges[i] == points[j];
            if (flag) edgesIndices.push_back(j);
        }
        if (!flag) {
			points.push_back(edges[i]);
            edgesIndices.push_back(points.size() - 1);
        }
    }
	this->points = points;
    std::vector<std::set<uint> > edgeList(points.size());
    for (uint i = 0; i < edgesIndices.size(); i += 2) {
        uint ind1 = std::min(edgesIndices[i], edgesIndices[i + 1]);
        uint ind2 = std::max(edgesIndices[i], edgesIndices[i + 1]);
        if (ind1 == ind2) continue;
        edgeList[ind1].insert(ind2);
    }
    triangles.clear();
    for (uint i = 0; i < points.size(); ++i) {
        for (auto j = edgeList[i].begin(); j != edgeList[i].end(); ++j) {
            for (auto h = j; h != edgeList[i].end(); ++h) {
				if (*h == *j) continue;
                if (edgeList[*j].find(*h) != edgeList[*j].end()) {
                    triangles.push_back(i);
                    triangles.push_back(*j);
                    triangles.push_back(*h);
                }
            }
        }
    }
}

std::vector<vec4>& Geometry::getPoints() {
    return points;
}

std::vector<uint>& Geometry::getTriangles() {
    return triangles;
}

std::vector<uint>& Geometry::getEdges() {
    //return edges;
}

const std::vector<vec4>& Geometry::getPoints() const {
    const std::vector<vec4> &rf = points;
    return rf;
}

const std::vector<uint>& Geometry::getTriangles() const {
    const std::vector<uint> &rf = triangles;
    return rf;
}

const std::vector<uint>& Geometry::getEdges() const {
    //const std::vector<uint> &rf = edges;
    //return rf;
}

std::vector<vec4> Geometry::getEdgesAsPoints() const {
    std::vector<vec4> result;
    //for (uint i = 0; i < edges.size(); ++i)
    //    result.push_back(points[edges[i]]);
    return result;
}


uint Geometry::getTrianglesNumber() const {
    return triangles.size() / 3;
}

uint Geometry::getEdgesNumber() const {
//    return edges.size() / 2;
}

uint Geometry::getPointsNumber() const {
    return points.size();
}

vec4 Geometry::getTrianglePlant(const uint index) const {
    vec3 a = points[triangles[3 * index + 0]].xyz();
    vec3 b = points[triangles[3 * index + 1]].xyz();
    vec3 c = points[triangles[3 * index + 2]].xyz();
    vec3 norm = normalize(cross(b - a, c - a));
    return vec4(norm, dot(norm, a));
}

void Geometry::removeNullEdges() {
    std::vector<uint> newEdges;
//    for (uint i = 0; i < edges.size(); i += 2) {
  //      if (edges[i] == edges[i + 1]) continue;
    //    newEdges.push_back(edges[i + 0]);
      //  newEdges.push_back(edges[i + 1]);
    //}
    //edges = newEdges;
}

void Geometry::mergeSimilarPoints() {
	std::cout << "Merge similar points begin: " << points.size() << ' ' << triangles.size() << std::endl;
	std::vector<vec4> goodPoints;
	std::vector<uint> IndexMap(points.size());
	for (uint i = 0; i < points.size(); ++i) {
        bool flag = false;
        uint goodInd;
        for (uint j = 0; j < goodPoints.size() && !flag; ++j) {
            goodInd = j;
            flag = points[i] == goodPoints[j];
        }
        if (flag) IndexMap[i] = goodInd;
        else {
			IndexMap[i] = goodPoints.size();
            goodPoints.push_back(points[i]);
        }
	}
	points = goodPoints;
    std::vector<uint> newTri;
	for (uint i = 0; i < triangles.size(); i += 3) {
        triangles[i + 0] = IndexMap[triangles[i + 0]];
        triangles[i + 1] = IndexMap[triangles[i + 1]];
        triangles[i + 2] = IndexMap[triangles[i + 2]];
        sort(triangles.begin() + i, triangles.begin() + i + 3);
        if (triangles[i] == triangles[i + 1] || triangles[i + 1] == triangles[i + 2]) continue;
        newTri.push_back(triangles[i + 0]);
        newTri.push_back(triangles[i + 1]);
        newTri.push_back(triangles[i + 2]);
	}
	triangles = newTri;
	std::cout << "Merge similar points end: " << points.size() << ' ' << triangles.size() << std::endl;
}

void Geometry::removeBadTriangles() {
    std::vector<uint> newTriangles;
    for (uint i = 0; i < triangles.size(); i += 3) {
        sort(triangles.begin() + i, triangles.begin() + i + 3);
        if (triangles[i] == triangles[i + 1] || triangles[i + 1] == triangles[i + 2]) continue;
        newTriangles.push_back(triangles[i + 0]);
        newTriangles.push_back(triangles[i + 1]);
        newTriangles.push_back(triangles[i + 2]);
    }
    triangles = newTriangles;
}


void Geometry::removeSmallTriangles(const float minSqure) {
	std::cout << "Remove small triangles begin: " << points.size() << ' ' << triangles.size() << std::endl;
	uint cnt = 0;
    for (uint i = 0; i < triangles.size(); i += 3) {
		if (100 * (i + 3) / triangles.size() > 100 * i / triangles.size())
			std::cout << 100 * i / triangles.size() << "% small triangles removed " << cnt << std::endl;
        vec4 a = points[triangles[i + 0]];
		vec4 b = points[triangles[i + 1]];
		vec4 c = points[triangles[i + 2]];
        if (length(cross((a - b).xyz(), (a - c).xyz()) / 2) < minSqure / 100) {
			cnt++;
            mergeTwoPoints(triangles[i], triangles[i + 1], (a + b + c) / 3);
            mergeTwoPoints(triangles[i], triangles[i + 2], (a + b + c) / 3);
        }
    }
    std::cout << "Remove small triangles end: " << points.size() << ' ' << triangles.size() << std::endl;
}


void Geometry::removeEqualEdges() {
    /*std::vector<std::set<uint> > edgeList(points.size());
	for (uint i = 0; i < edges.size(); i += 2) {
        sort(edges.begin() + i, edges.begin() + i + 2);
		edgeList[edges[i]].insert(edges[i + 1]);
	}
	std::vector<uint> newEdges;
	for (uint i = 0; i < edgeList.size(); ++i) {
        for (auto j = edgeList[i].begin(); j != edgeList[i].end(); ++j) {
            newEdges.push_back(i);
            newEdges.push_back(*j);
        }
	}
	edges = newEdges;*/
}

void Geometry::repairGeometry() {
	/*std::cout << points.size() << ' ' << edges.size() << ' ' << triangles.size() << std::endl;
    mergeSimilarPoints();
    std::cout << points.size() << ' ' << edges.size() << ' ' << triangles.size() << std::endl;
    removeNullEdges();
    std::cout << points.size() << ' ' << edges.size() << ' ' << triangles.size() << std::endl;
    removeEqualEdges();
    std::cout << points.size() << ' ' << edges.size() << ' ' << triangles.size() << std::endl;
    removeBadTriangles();
    std::cout << points.size() << ' ' << edges.size() << ' ' << triangles.size() << std::endl;*/
}


void Geometry::Write(const std::string &path) const {
    std::ofstream out(path);
	out << points.size() << std::endl;
	out.precision(10);
    for (auto v: points)
        out << std::fixed << v.x << ' ' << std::fixed << v.y << ' ' <<
			   std::fixed << v.z << ' ' << std::fixed << v.w << ' ';
	out << std::endl << triangles.size() << std::endl;
    for (auto i: triangles)
        out << i << ' ';
    out.close();
}

void Geometry::Read(const std::string &path) {
    std::ifstream in(path);
    uint pC;
    in >> pC;
    points.resize(pC);
    for (uint i = 0; i < pC; ++i) {
        vec4 v;
        in >> v.x >> v.y >> v.z >> v.w;
        points[i] = v;
    }
    in >> pC;
    triangles.resize(pC);
    for (uint i = 0; i < pC; ++i)
        in >> triangles[i];
	in.close();
}

std::vector<std::set<uint> > Geometry::getAdjacencyListIndices() {
	std::vector<std::set<uint> > result(points.size());
    for (uint i = 0; i < triangles.size(); i += 3) {
        std::sort(triangles.begin() + i, triangles.begin() + i + 3);
        result[triangles[i + 0]].insert(triangles[i + 1]);
        result[triangles[i + 0]].insert(triangles[i + 2]);
        result[triangles[i + 1]].insert(triangles[i + 2]);
    }
    return result;
}

std::vector<uint> Geometry::getEdgeListIndices(){
	std::vector<uint> result;
    auto adjListInd = getAdjacencyListIndices();
    for (uint i = 0; i < adjListInd.size(); ++i) {
        for (auto j = adjListInd[i].begin(); j != adjListInd[i].end(); ++j) {
            result.push_back(i);
            result.push_back(*j);
        }
    }
    return result;
}

std::vector<vec4> Geometry::getEdgeList(){
	std::vector<vec4> result;
	auto edgeList = getEdgeListIndices();
    for (uint i = 0; i < edgeList.size(); ++i) {
        result.push_back(points[edgeList[i]]);
    }
    return result;
}

std::vector<std::vector<vec4> > Geometry::getOutLines() {
    std::vector<std::vector<vec4> > result(points.size());
    auto edgeList = getAdjacencyListIndices();
    for (uint i = 0; i < points.size(); ++i) {
        for (auto j = edgeList[i].begin(); j != edgeList[i].end(); ++j) {
            vec4 newVec = points[*j] - points[i];
			bool flag = true;
            for (uint h = 0; h < result[i].size() && flag; ++h)
                flag = length(cross(newVec.xyz(), result[i][h].xyz())) >= VEC_EPS;
			if (flag)
                result[i].push_back(newVec);
        }
    }
    return result;
}

std::vector<std::set<uint> > Geometry::getBidirAdjacencyListIndices() {
	std::vector<std::set<uint> > result(points.size());
    for (uint i = 0; i < triangles.size(); i += 3) {
        result[triangles[i + 0]].insert(triangles[i + 1]);
        result[triangles[i + 0]].insert(triangles[i + 2]);
        result[triangles[i + 1]].insert(triangles[i + 0]);
        result[triangles[i + 1]].insert(triangles[i + 2]);
        result[triangles[i + 2]].insert(triangles[i + 0]);
        result[triangles[i + 2]].insert(triangles[i + 1]);
    }
    return result;
}

std::vector<uint> splitPolygonByTriangles(const std::vector<uint> &polygon) {
    //std::cout << "Split polygons by triangles begin" << std::endl;
    std::vector<uint> result;
    for (uint i = 2; i < polygon.size(); ++i) {
        result.push_back(polygon[0]);
        result.push_back(polygon[i - 1]);
        result.push_back(polygon[i]);
    }
    //std::cout << "Split polygons by triangles end" << std::endl;
    return result;
}


std::vector<std::vector<uint> > mergeTriangles(const std::vector<std::vector<uint> > &triangles) {
    //std::cout << "Merge triangles begin" << std::endl;
    std::map<uint, std::map<uint, uint> > neibors;
    for (uint i = 0; i < triangles.size(); i++) {
        neibors[triangles[i][0]][triangles[i][1]]++;
        neibors[triangles[i][0]][triangles[i][2]]++;
        neibors[triangles[i][1]][triangles[i][0]]++;
        neibors[triangles[i][1]][triangles[i][2]]++;
        neibors[triangles[i][2]][triangles[i][1]]++;
        neibors[triangles[i][2]][triangles[i][0]]++;
    }
    std::set<uint> cont;
    for (auto i = neibors.begin(); i != neibors.end(); ++i) {
        uint sum = 0;
        for (auto j = i->second.begin(); j != i->second.end(); ++j)
            sum += j->second & 1;
        if (sum) cont.insert(i->first);
    }
    std::map<uint, std::set<uint> > edgeList;
    for (uint i = 0; i < triangles.size(); ++i) {
        std::vector<uint> subcont;
        for (uint j = 0; j < 3; ++j)
            if (cont.find(triangles[i][j]) != cont.end())
                subcont.push_back(triangles[i][j]);
        for (uint j = 0; j < subcont.size(); ++j) {
            for (uint h = 0; h < subcont.size(); ++h) {
                if (h == j) continue;
                edgeList[triangles[i][j]].insert(triangles[i][h]);
            }
        }
    }
    std::vector<std::vector<uint> > result;
    std::set<uint> removed;
    for (auto i = cont.begin(); i != cont.end(); ++i) {
        uint value = *i;
        if (removed.find(value) != removed.end()) continue;
        result.push_back(std::vector<uint>());
        while (removed.find(value) == removed.end()) {
            result.back().push_back(value);
            removed.insert(value);
            for (auto j = edgeList[value].begin(); j != edgeList[value].end(); ++j)
                if (removed.find(*j) == removed.end() && cont.find(*j) != cont.end()) {
                    value = *j;
                    break;
                }
        }
    }
    //std::cout << "Merge triangles end" << std::endl;
    return result;
}

void Geometry::removeExcessPoints() {
    std::cout << "Remove excess points begin" << std::endl;
    std::set<uint> badPoints;
    for (uint i = 0; i < points.size(); ++i)
        badPoints.insert(i);
    for (uint i = 0; i < triangles.size(); ++i)
        badPoints.erase(triangles[i]);
    std::vector<vec4> newPoints;
    std::vector<uint> newInd(points.size());
    for (uint i = 0; i < points.size(); ++i) {
        if (badPoints.find(i) == badPoints.end()) {
            newInd[i] = newPoints.size();
            newPoints.push_back(points[i]);
        }
    }
    points = newPoints;
    for (uint i = 0; i < triangles.size(); ++i)
        triangles[i] = newInd[triangles[i]];
    std::cout << "Remove excess points end" << std::endl;
}

std::vector<std::vector<uint> > removeNotAnglePoints(const std::vector<std::vector<uint> > &polys, const std::vector<vec4> points) {
    std::vector<std::vector<uint> > res(polys.size());
    for (uint i = 0; i < polys.size(); ++i) {
        vec3 direct = (points[polys[i][0]] - points[polys[i][1]]).xyz();
        for (uint j = 0; j < polys[i].size(); ++j) {
            vec3 newDir = (points[polys[i][(j + 1) % polys[i].size()]] - points[polys[i][j]]).xyz();
            if (length(cross(direct, newDir)) > VEC_EPS) {
                res[i].push_back(polys[i][j]);
            }
            direct = newDir;
        }
    }
    return res;
}

std::vector<std::vector<std::vector<uint> > > Geometry::groupTrianglesByPlanes() {
    std::cout << "Group triangles by planes begin" << std::endl;
    std::vector<std::vector<std::vector<uint> > > result;
    std::vector<vec4> planes;
    for (uint i = 0; i < triangles.size(); i += 3) {
        vec4 plane = getTrianglePlant(i / 3);
        uint ind = 0;
        for (; ind < planes.size() && planes[ind] != plane; ++ind);
        if (ind == planes.size()) {
            planes.push_back(plane);
            result.resize(result.size() + 1);
        }
        result[ind].push_back(std::vector<uint>(triangles.begin() + i, triangles.begin() + i + 3));
    }
    std::cout << "Group triangles by planes end" << std::endl;
    return result;
}

void Geometry::removeInternalPoints() {
    std::cout << "Remove internal points begin" << std::endl;
    std::vector<uint> newTriangles;
    for (auto &plant: groupTrianglesByPlanes()) {
        for (auto &polygon: mergeTriangles(plant)) {
            for (auto i: splitPolygonByTriangles(polygon)) {
                newTriangles.push_back(i);
            }
        }
        break;
    }
    triangles = newTriangles;
    removeExcessPoints();
    std::cout << "Remove internal points end" << std::endl;
}

void Geometry::FirstPlane() {
    auto plant = groupTrianglesByPlanes()[0];
    std::vector<uint> newTriangles;
    for (uint i = 0; i < plant.size(); ++i)
        for (uint j = 0; j < plant[i].size(); ++j)
            newTriangles.push_back(plant[i][j]);
    triangles = newTriangles;
    removeExcessPoints();
}
