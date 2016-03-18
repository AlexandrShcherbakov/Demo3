#include "Simplifier.h"
#include "Geometry.h"

///Convenient small functions
geometry make_geometry(const std::vector<vec3> points, const std::vector<uint> indices) {
	return std::make_tuple(points, indices);
}

///*****************************************************************************
///QuadricErrorMetricsSimplifier
///*****************************************************************************

/*geometry prepareGeometry(const geometry baseGeometry) {
	auto indices = std::get<1>(baseGeometry);
	auto points = std::get<0>(baseGeometry);

	std::set<vec3> newPointsSet;
	for (auto &p: points)
		newPointsSet.insert(p);

	std::map<vec3, uint> pointsIndides;
	std::vector<vec3> pointsNew;
	for (auto &p: newPointsSet) {
		pointsIndides[p] = pointsNew.size();
		pointsNew.push_back(p);
	}

	std::vector<uint> triangles;
	for (uint i = 0; i < points.size(); i += 3) {
		vec3 a = points[i + 0];
		vec3 b = points[i + 1];
		vec3 c = points[i + 2];
		if (length(cross(b - a, c - a)) == 0)
			continue;
		triangles.push_back(pointsIndides[a]);
		triangles.push_back(pointsIndides[b]);
		triangles.push_back(pointsIndides[c]);
	}

	return make_tuple(pointsNew, triangles);
}

geometry removeInternalPoints(const geometry baseGeometry) {
	std::vector<uint> indices;
	for (uint i = 0; i < std::get<1>(baseGeometry).size(); ++i)
		indices.push_back(std::get<1>(baseGeometry)[i]);
	std::vector<vec3> points;
	for (uint i = 0; i < std::get<0>(baseGeometry).size(); ++i)
		points.push_back(std::get<0>(baseGeometry)[i]);

	std::map<vec3, std::set<vec3> > normals;
	for (uint i = 0; i < indices.size(); i += 3) {
		vec3 a = points[indices[i + 0]];
		vec3 b = points[indices[i + 1]];
		vec3 c = points[indices[i + 2]];
		if (length(cross(b - a, c - a)) == 0)
			continue;
		vec3 norm = normalize(cross(b - a, c - a));
		normals[a].insert(norm);
		normals[b].insert(norm);
		normals[c].insert(norm);
	}

	for (int i = points.size() - 1; i >= 0; --i) {
		if (100 * (points.size() - i) / points.size() > 100 * (points.size() - i - 1) / points.size())
			std::cout << 100 * (points.size() - i) / points.size() << "% remove internal points process" << std::endl;
		if (normals[points[i]].size() > 1) continue;
		vec3 oldPoint = points[i];
		int j;
		for (j = 0; j < indices.size() && indices[j] != i; ++j);
		vec3 newTarget;
		uint tarNum;
		if (j % 3 != 2) {
			newTarget = points[indices[j + 1]];
			tarNum = indices[j + 1];
		}
		else {
			newTarget = points[indices[j - 1]];
			tarNum = indices[j - 1];
		}
		for (j = indices.size() - 3; j > 0; j -= 3) {
			vec3 a = points[indices[j + 0]];
			vec3 b = points[indices[j + 1]];
			vec3 c = points[indices[j + 2]];
			if (oldPoint == a || oldPoint == b || oldPoint == c) {
				if (newTarget == a || newTarget == b || newTarget == c) {
					indices.erase(indices.begin() + j, indices.begin() + j + 3);
				} else {
					for (uint h = 0; h < 3; ++h)
						if (points[indices[j + h]] == oldPoint)
							indices[j + h] = tarNum;
				}
			}
		}
	}
	return prepareGeometry(make_tuple(points, indices));
}


std::vector<mat4> computePlantMatrices(const geometry baseGeometry) {
	auto indices = std::get<1>(baseGeometry);
	auto points = std::get<0>(baseGeometry);
	std::vector<mat4> result(indices.size() / 3);
	for (uint i = 0; i < indices.size(); i += 3) {
		vec3 a = points[indices[i + 0]];
		vec3 b = points[indices[i + 1]];
		vec3 c = points[indices[i + 2]];

		vec3 norm = normalize(cross(b - a, c - a));
		float d = -dot(norm, a);
		vec4 plant(norm, d);
		mat4 K;
		for (uint j = 0; j < 4u; ++j)
			for (uint h = 0; h < 4u; ++h)
				K[j][h] = plant[j] * plant[h];
		result[i / 3] = K;
	}
	return result;
}

std::vector<mat4> computePointMatrices(const geometry baseGeometry, const std::vector<mat4> matrices) {
	auto indices = std::get<1>(baseGeometry);
	auto points = std::get<0>(baseGeometry);
	std::vector<mat4> result(points.size());
	for (uint i = 0; i < indices.size(); i++) {
		result[indices[i]] += matrices[i / 3];
	}
	return result;
}*/

class QEMEdge {
private:
	void buildTarget() {
		targetmat = startmat[0] + startmat[1];
		/*mat4 tmpmat = targetmat;
		tmpmat[3][0] = tmpmat[3][1] = tmpmat[3][2] = 0.0f;
		tmpmat[3][3] = 1.0f;
		targetvec = tmpmat.unmatrixN3() * vec4(0, 0, 0, 1);*/
		targetvec = (startvec[0] + startvec[1]) / 2;
		//error = dot(targetmat * targetvec, targetvec);
		error = length(startvec[0] - startvec[1]);
	}
public:
	vec4 startvec[2], targetvec;
	mat4 startmat[2], targetmat;
	float error;

	QEMEdge() {};

	QEMEdge(vec4 v1, vec4 v2, mat4 m1, mat4 m2) {
		startvec[0] = v1;
		startvec[1] = v2;
		startmat[0] = m1;
		startmat[1] = m2;
		buildTarget();
	}

	bool operator < (const QEMEdge& e1) const {
		return error < e1.error;
	}

	bool operator > (const QEMEdge& e1) const {
		return error > e1.error;
	}

	bool operator == (const QEMEdge& e1) const {
		return error == e1.error;
	}

	void update(const QEMEdge& e) {
	    bool flag = false;
		for (uint i = 0; i < 2; ++i) {
			if (startvec[i] == e.startvec[0] || startvec[i] == e.startvec[1]) {
				startvec[i] = e.targetvec;
				startmat[i] = e.targetmat;
				flag = true;
			}
		}
		if (flag)
            buildTarget();
	}
};

/*void heap_push(std::vector<QEMEdge> &heap, QEMEdge x) {
	heap.push_back(x);
	for (uint i = heap.size() - 1; i > 1 && heap[i] < heap[i / 2]; i /= 2)
        std::swap(heap[i], heap[i / 2]);
}


std::vector<QEMEdge> createHeapOfEdges(const geometry baseGeometry, const std::vector<mat4> pointMat) {
	auto indices = std::get<1>(baseGeometry);
	auto points = std::get<0>(baseGeometry);
	std::map<uint, std::set<uint> > edges;
	for (uint i = 0; i < indices.size(); i += 3) {
		sort(indices.begin() + i, indices.begin() + i + 3);
		edges[indices[i + 0]].insert(indices[i + 1]);
		edges[indices[i + 0]].insert(indices[i + 2]);
		edges[indices[i + 1]].insert(indices[i + 2]);
	}

    std::vector<QEMEdge> edgesHeap(1);
    for (auto i = edges.begin(); i != edges.end(); ++i) {
		for (auto j = i->second.begin(); j != i->second.end(); ++j) {
			vec4 v1 = vec4(points[i->first], 1);
			vec4 v2 = vec4(points[*j], 1);
			mat4 m1 = pointMat[i->first];
			mat4 m2 = pointMat[*j];
			heap_push(edgesHeap, QEMEdge(v1, v2, m1, m2));
		}
    }
	return edgesHeap;
}

void heapify(std::vector<QEMEdge>& heap) {
    for (uint i = 1; i < heap.size() / 2;) {
        if (heap[i] > heap[i * 2]) {
            if (i * 2 + 1 < heap.size() && heap[i * 2] > heap[i * 2 + 1]) {
				std::swap(heap[i], heap[i * 2 + 1]);
				i = 2 * i + 1;
            } else {
				std::swap(heap[i], heap[i * 2]);
				i = 2 * i;
            }
        } else if (i * 2 + 1 < heap.size() && heap[i] > heap[i * 2 + 1]) {
        	std::swap(heap[i], heap[i * 2 + 1]);
			i = 2 * i + 1;
        } else {
            break;
        }
    }
}

void heap_updateElem(std::vector<QEMEdge>& heap, uint index) {
    if (index > 1 && heap[index / 2] > heap[index]) {
        for (; index > 1 && heap[index / 2] > heap[index]; index /= 2)
            std::swap(heap[index / 2], heap[index]);
    } else {
        for (; index < heap.size() / 2;) {
            float errP = heap[index].error;
            float errL = heap[index * 2].error;
            float errR = (index * 2 + 1) < heap.size() ? heap[index * 2 + 1].error : FLT_MAX;
            if (errP > errL) {
                if (errL > errR) {
                    std::swap(heap[index], heap[index * 2 + 1]);
                    index = 2 * index + 1;
                } else {
                    std::swap(heap[index], heap[index * 2]);
                    index = 2 * index;
                }
            } else if (errP > errR) {
                std::swap(heap[index], heap[index * 2 + 1]);
                index = 2 * index + 1;
            } else {
                break;
            }
        }
    }
}

void removeSmallEdges(std::vector<QEMEdge>& heap, const uint removeIters) {
	std::cout << heap.size() << ' ' << removeIters << std::endl;
	for (uint i = 0; i < removeIters; ++i) {
		if (100 * (i + 1) / removeIters > 100 * (i) / removeIters)
			std::cout << 100 * (i + 1) / removeIters << "% remove edges process" << std::endl;
		QEMEdge updated = heap[1];
		std::swap(heap[1], heap[heap.size() - 1]);
		heap.pop_back();
		heapify(heap);
		for (uint j = 1; j < heap.size(); ++j) {
			heap[i].update(updated);
			heap_updateElem(heap, i);
		}
	}
}

geometry createNewGeometry(std::vector<QEMEdge>& heap) {
    std::cout << heap.size() << std::endl;
	std::set<vec4> pointsSet;
	for (auto &edge: heap) {
		pointsSet.insert(edge.startvec[0]);
		pointsSet.insert(edge.startvec[1]);
	}
	std::vector<vec3> pointsVec;
	std::map<vec4, uint> pointsInd;
	for (auto &point: pointsSet) {
		pointsVec.push_back(vec3(point.x, point.y, point.z));
		pointsInd[point] = pointsVec.size() - 1;
	}

	std::vector<std::set<uint> > edges(pointsInd.size());
	for (auto &edge: heap) {
		uint i1 = std::min(pointsInd[edge.startvec[0]], pointsInd[edge.startvec[1]]);
		uint i2 = std::max(pointsInd[edge.startvec[0]], pointsInd[edge.startvec[1]]);
		edges[i1].insert(i2);
	}

	std::vector<uint> triangles;
	for (uint i = 0; i < edges.size(); ++i) {
		for (auto j = edges[i].begin(); j != edges[i].end(); ++j) {
			std::vector<uint> third(edges[i].size() + edges[*j].size());
			auto end = std::set_intersection(edges[i].begin(), edges[i].end(), edges[*j].begin(), edges[*j].end(), third.begin());
			third.resize(end - third.begin());
			for (auto &h : third) {
				triangles.push_back(i);
				triangles.push_back(*j);
				triangles.push_back(h);
			}
		}
	}
	std::cout << pointsVec.size() << ' ' << triangles.size() << std::endl;
	return std::make_tuple(pointsVec, triangles);
}

geometry QuadricErrorMetricsSimplifier(const geometry& inGeometry, const uint polygonLimit) {
	///Step 0: prepare data
	auto baseGeometry = prepareGeometry(inGeometry); std::cout << "Step 0 complete\n";
	///Step 1: remove internal points
	baseGeometry = removeInternalPoints(baseGeometry); std::cout << "Step 1 complete\n";
	///Step 2: compute error-matrices for plants
	auto plantMat = computePlantMatrices(baseGeometry); std::cout << "Step 2 complete\n";
	///Step 3: compute error-matrices for points
	auto pointMat = computePointMatrices(baseGeometry, plantMat); std::cout << "Step 3 complete\n";
	///Step 4: create heap with edges
	auto heap = createHeapOfEdges(baseGeometry, pointMat); std::cout << "Step 4 complete\n";
	///Step 5: remove some edges
	removeSmallEdges(heap, (std::get<1>(baseGeometry).size() - polygonLimit)); std::cout << "Step 5 complete\n";
	///Step 6: combine new scene from leftover points
	return createNewGeometry(heap);
}*/

std::pair<std::vector<vec4>, std::vector<mat4> > NewcomputePlantMatrices(const Geometry &geom) {
    std::map<vec4, mat4> result;
    std::vector<vec4> resVecs;
    std::vector<mat4> resMats;
    std::cout << geom.getTrianglesNumber() << std::endl;
	for (uint i = 0; i < geom.getTrianglesNumber(); ++i) {
        vec4 plant = geom.getTrianglePlant(i);
        bool flag = false;
        for (uint j = 0; j < resVecs.size() && !flag; ++j) {
            flag = flag || plant == resVecs[j];
        }
        if (flag) continue;
		mat4 K;
		for (uint j = 0; j < 4u; ++j)
			for (uint h = 0; h < 4u; ++h)
				K[j][h] = plant[j] * plant[h];
        resVecs.push_back(plant);
        resMats.push_back(K);
	}
	return std::make_pair(resVecs, resMats);
}

std::vector<mat4> NewcomputePointMatrices(const Geometry& geom, std::pair<std::vector<vec4>, std::vector<mat4> > matrices) {
    std::vector<std::vector<vec4> > plantsForPoints(geom.getPointsNumber());
    for (uint i = 0; i < geom.getTrianglesNumber(); ++i) {
        uint a = geom.getTriangles()[3 * i + 0];
        uint b = geom.getTriangles()[3 * i + 1];
        uint c = geom.getTriangles()[3 * i + 2];
        vec4 plant = geom.getTrianglePlant(i);
        bool flag = false;
        for (uint j = 0; j < plantsForPoints[a].size() && !flag; ++j) {
            flag = flag || plantsForPoints[a][j] == plant;
        }
        if (!flag) plantsForPoints[a].push_back(plant);
        flag = false;
        for (uint j = 0; j < plantsForPoints[b].size() && !flag; ++j) {
            flag = flag || plantsForPoints[b][j] == plant;
        }
        if (!flag) plantsForPoints[b].push_back(plant);
        flag = false;
        for (uint j = 0; j < plantsForPoints[c].size() && !flag; ++j) {
            flag = flag || plantsForPoints[c][j] == plant;
        }
        if (!flag) plantsForPoints[c].push_back(plant);
    }
    std::vector<mat4> result(geom.getPointsNumber());
    for (uint i = 0; i < geom.getPointsNumber(); ++i) {
        for (auto j = 0; j != plantsForPoints[i].size(); ++j) {
            for (uint h = 0; h < matrices.first.size(); ++i) {
                if (matrices.first[h] == plantsForPoints[i][j]) {
                    result[i] += matrices.second[h];
                }
            }
        }
    }
	return result;
}

void heapify(std::vector<QEMEdge>& heap, uint i) {
    while (2 * i < heap.size()) {
		if (heap[i] > heap[2 * i]) {
			if (2 * i + 1 < heap.size() && heap[2 * i] > heap[2 * i + 1]) {
				std::swap(heap[i], heap[2 * i + 1]);
				i = 2 * i + 1;
			} else {
                std::swap(heap[i], heap[2 * i]);
                i = 2 * i;
			}
		} else {
			if (2 * i + 1 < heap.size() && heap[i] > heap[2 * i + 1]) {
				std::swap(heap[i], heap[2 * i + 1]);
				i = 2 * i + 1;
			} else {
                break;
			}
		}
    }
}

std::vector<QEMEdge> NewcreateHeapOfEdges(const Geometry& geom, const std::vector<mat4> pointMat) {
	std::vector<QEMEdge> edgesHeap;
	for (uint i = 0; i < geom.getEdgesNumber(); ++i) {
        if (100 * (i + 1) / geom.getEdgesNumber() > 100 * i / geom.getEdgesNumber())
            std::cout << 100 * i / geom.getEdgesNumber() << std::endl;
        vec4 v1 = geom.getPoints()[geom.getEdges()[2 * i + 0]];
        vec4 v2 = geom.getPoints()[geom.getEdges()[2 * i + 1]];
        mat4 m1 = pointMat[geom.getEdges()[2 * i + 0]];
        mat4 m2 = pointMat[geom.getEdges()[2 * i + 1]];
        edgesHeap.push_back(QEMEdge(v1, v2, m1, m2));
	}
	std::cout << "Start sort" << std::endl;
	sort(edgesHeap.rbegin(), edgesHeap.rend());
	//sort(edgesHeap.begin(), edgesHeap.end());
	return edgesHeap;
}

void NewremoveSmallEdges(std::vector<QEMEdge>& heap, const uint removeIters) {
	std::cout << heap.size() << ' ' << removeIters << std::endl;
    std::list<QEMEdge> heapList(heap.begin(), heap.end());
	for (uint i = 0; i < removeIters; ++i) {
		if (100 * (i + 1) / removeIters > 100 * (i) / removeIters)
			std::cout << 100 * (i + 1) / removeIters << "% remove edges process " << std::endl;
		QEMEdge updated = heap.back();
		//std::cout << updated.startvec[0] << ' ' << updated.startvec[1] << ' ' << updated.targetvec << std::endl;
		heap.pop_back();
		uint swaps = 0;
		std::vector<QEMEdge> updList;
		//std::cout << updated.error << std::endl;
        for (auto j = heapList.begin(); j != heapList.end(); ++j) {
            float old = j->error;
			j->update(updated);
            if (j->error != old) {
                //std::cout << old << ' ' << j->error << ' ' << std::endl;
                updList.push_back(*j);
                j = heapList.erase(j);
                j--;
                /*if (old > heap[j].error) {
                    for (uint h = j; h + 1 < heap.size() && heap[h] < heap[h + 1]; ++h) {
                        std::swap(heap[h], heap[h + 1]);
                        swaps++;
                    }
                } else {
                    for (uint h = j; h - 1 >= 0 && heap[h] > heap[h - 1]; --h) {
                        std::swap(heap[h], heap[h - 1]);
                        swaps++;
                    }
                }*/
            }
		}
		std::sort(updList.rbegin(), updList.rend());
        for (uint j = 1; j < updList.size(); ++j) {
            if (updList[j] > updList[j - 1]) std::cout << i << " bad updList order" << std::endl;
        }
		//std::sort(updList.begin(), updList.end());
		auto it1 = heapList.begin();
        auto it2 = it1++;
        uint j = 0;
        while (j < updList.size() && updList[j] > *it2) {
            heapList.insert(it2, updList[j]);
            j++;
        }
        while (j < updList.size() && it1 != heapList.end()) {
            if (*it2 > updList[j] && *it1 < updList[j] || *it2 == updList[j]) {
                it2 = heapList.insert(it1, updList[j]);
                //it2++;
                j++;
            } else {
                it1++;
                it2++;
            }
        }
        /*std::cout << j << ' ' << updList.size() << std::endl;
        it1 = heapList.end();
        for (uint j = 0; j < 10; ++j) it1--;
        while (it1 != heapList.end()) {
			std::cout << it1->error << std::endl;
			it1++;
        }
        for (auto x: updList) std::cout << x.error << ' ';
        std::cout << std::endl;
        std::cout << j << ' ' << updList.size() << std::endl;*/
        while (j < updList.size()) {
            heapList.push_back(updList[j++]);
        }
        /*for (uint j = 0; j < 10; ++j) it1--;
        while (it1 != heapList.end()) {
			std::cout << it1->error << std::endl;
			it1++;
        }
        std::cin >>j;*/
        //std::list<QEMEdge> updToList(updList.begin(), updList.end());
		//std::list<QEMEdge> newList(heapList.size() + updList.size());
        //std::merge(heapList.rbegin(), heapList.rend(), updToList.rbegin(), updToList.rend(), newList.begin());
        //heapList = newList;
		//std::cout << swaps << std::endl;
        //sort(heap.rbegin(), heap.rend());
	}
}

Geometry NewcreateNewGeometry(std::vector<QEMEdge>& heap) {
    std::cout << heap.size() << std::endl;
	Geometry geom;
	std::vector<vec4> edges;
	for (uint i = 0; i < heap.size(); ++i) {
        edges.push_back(heap[i].startvec[0]);
        edges.push_back(heap[i].startvec[1]);
	}
	geom.loadFromEdges(edges);
	return geom;
}

mat4 createMat4FromVec4(const vec4& v) {
    mat4 res;
    for (uint i = 0; i < 4; ++i)
        for (uint j = 0; j < 4; ++j)
            res[i][j] = v[i] * v[j];
    return res;
}

std::vector<mat4> computeMatrices(const Geometry& geom) {
    std::vector<mat4> result(geom.getPointsNumber());
    std::vector<std::vector<vec4> > plants(geom.getPointsNumber());
    for (uint i = 0; i < geom.getTrianglesNumber(); ++i) {
        if (100 * (i + 1) / geom.getTrianglesNumber() > 100 * i / geom.getTrianglesNumber())
            std::cout << 100 * i / geom.getTrianglesNumber() << "% triangles worked" << std::endl;
        vec4 plant = geom.getTrianglePlant(i);
        for (uint j = i * 3; j < 3 * i + 3; ++j) {
            uint index = geom.getTriangles()[j];
            bool flag = false;
            for (uint h = 0; h < plants[index].size() && !flag; ++h) {
                flag = flag || plants[index][h] == plant;
            }
            if (!flag) {
                result[index] += createMat4FromVec4(plant);
                plants[index].push_back(plant);
            }
        }
    }
    return result;
}

void WriteGeometry1(const Geometry & geom, const std::string &path) {
    std::ofstream out(path);
    out << geom.getPointsNumber() << std::endl;
    for (uint i = 0; i < geom.getPointsNumber(); ++i)
        for (uint j = 0; j < 4; ++j)
            out << geom.getPoints()[i][j] << ' ';
    out << geom.getTrianglesNumber() << std::endl;
    for (uint i = 0; i < geom.getTrianglesNumber() * 3; ++i)
        out << geom.getTriangles()[i] << ' ';
    out.close();
    std::cout << "OK" << std::endl;
}

Geometry NewQuadricErrorMetricsSimplifier(const Geometry inGeometry, const uint polygonLimit) {
    Geometry baseGeometry = inGeometry;
    ///Step 1: remove internal points
    //baseGeometry.removePointsIntoLines(); std::cout << "Step 0 complete\n";
    //baseGeometry.removePointsIntoPolygons(); std::cout << "Step 1 complete\n";
    //baseGeometry.removeSmallTriangles(VEC_EPS); std::cout << "Step 2 complete\n";
    //baseGeometry.Write("Before QEM");
	///Step 2: compute error-matrices for plants
	//auto plantMat = NewcomputePlantMatrices(baseGeometry); std::cout << "Step 2 complete\n";
	///Step 3: compute error-matrices for points
	//auto pointMat = NewcomputePointMatrices(baseGeometry, plantMat); std::cout << "Step 3 complete\n";
    auto pointMat = computeMatrices(baseGeometry); std::cout << "Step 3 complete" << std::endl;
	///Step 4: create heap with edges
	auto heap = NewcreateHeapOfEdges(baseGeometry, pointMat); std::cout << "Step 4 complete\n";
	///Step 5: remove some edges
	NewremoveSmallEdges(heap, (baseGeometry.getTrianglesNumber() - polygonLimit) * 3 / 2); std::cout << "Step 5 complete\n";
	//NewremoveSmallEdges(heap, polygonLimit); std::cout << "Step 5 complete\n";
	//NewremoveSmallEdges(heap, 1); std::cout << "Step 5 complete\n";
	///Step 6: combine new scene from leftover points
	return NewcreateNewGeometry(heap);
}

///*****************************************************************************
///SmallEdgesSimplifier
///*****************************************************************************


Geometry SmallEdgesSimplifier(const Geometry &inGeometry, const uint polygonLimit) {
    Geometry baseGeometry = inGeometry;
    //std::cout << baseGeometry.getPointsNumber() << ' ' << baseGeometry.getEdgesNumber() << ' ' << baseGeometry.getTrianglesNumber() << std::endl;
    for (uint i = 0; i < polygonLimit; ++i) {
        if (100 * (i + 1) / polygonLimit > 100 * i / polygonLimit)
			std::cout << 100 * i / polygonLimit << "% edges removed" << std::endl;
		uint mn = 0;
		vec4 target;
		float err = 1e9;
        auto edgeList = baseGeometry.getEdgeList();
        for (uint j = 0; j < edgeList.size(); j += 2) {
            vec4 a = edgeList[j];
            vec4 b = edgeList[j + 1];
            if (length(a - b) < err) {
                err = length(a - b);
                mn = j;
                target = (a + b) / 2;
            }
        }
        baseGeometry.mergeTwoPoints(baseGeometry.getEdges()[mn + 0],
									baseGeometry.getEdges()[mn + 1],
									target);
    }
    return baseGeometry;
}


///*****************************************************************************
///SOMmapSimplifier
///*****************************************************************************


geometry SOMmapSimplifier(const geometry& baseGeometry, const uint polygonLimit) {
	///Step 1: create initial network
	///Step 2: prepare functions for training
	///Step 3: train network
	///Step 4: combine new scene from network
}
