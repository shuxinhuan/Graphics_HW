#include "mesh.h"
#include <iostream>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <queue>


HEdge::HEdge(bool b) {
	mBoundary = b;

	mTwin = nullptr;
	mPrev = nullptr;
	mNext = nullptr;

	mStart = nullptr;
	mFace = nullptr;

	mFlag = false;
	mValid = true;
}

HEdge* HEdge::twin() const {
	return mTwin;
}

HEdge* HEdge::setTwin(HEdge* e) {
	mTwin = e;
	return mTwin;
}

HEdge* HEdge::prev() const {
	return mPrev;
}

HEdge* HEdge::setPrev(HEdge* e) {
	mPrev = e;
	return mPrev;
}

HEdge* HEdge::next() const {
	return mNext;
}

HEdge* HEdge::setNext(HEdge* e) {
	mNext = e;
	return mNext;
}

Vertex* HEdge::start() const {
	return mStart;
}

Vertex* HEdge::setStart(Vertex* v) {
	mStart = v;
	return mStart;
}

Vertex* HEdge::end() const {
	return mNext->start();
}

Face* HEdge::leftFace() const {
	return mFace;
}

Face* HEdge::setFace(Face* f) {
	mFace = f;
	return mFace;
}

bool HEdge::flag() const {
	return mFlag;
}

bool HEdge::setFlag(bool b) {
	mFlag = b;
	return mFlag;
}

bool HEdge::isBoundary() const {
	return mBoundary;
}

bool HEdge::isValid() const {
	return mValid;
}

bool HEdge::setValid(bool b) {
	mValid = b;
	return mValid;
}

OneRingHEdge::OneRingHEdge(const Vertex* v) {
	if (v == nullptr) {
		mStart = nullptr;
		mNext = nullptr;
	} else {
		mStart = v->halfEdge();
		mNext = v->halfEdge();
	}
}

HEdge* OneRingHEdge::nextHEdge() {
	HEdge* ret = mNext;
	if (mNext != nullptr && mNext->prev()->twin() != mStart) {
		mNext = mNext->prev()->twin();
	} else {
		mNext = nullptr;
	}
	return ret;
}

OneRingVertex::OneRingVertex(const Vertex* v): ring(v) {
}

Vertex* OneRingVertex::nextVertex() {
	HEdge* he = ring.nextHEdge();
	return he != nullptr ? he->end() : nullptr;
}

Vertex::Vertex() : mHEdge(nullptr), mFlag(0), mKnown(false) {
	mPosition = Eigen::Vector3f::Zero();
	mColor = VCOLOR_BLUE;
	mNormal = Eigen::Vector3f::Zero();
}

Vertex::Vertex(const Eigen::Vector3f& v): mPosition(v), mHEdge(nullptr), mFlag(0), mKnown(false) {
	mColor = VCOLOR_BLUE;
	mNormal = Eigen::Vector3f::Zero();
}

Vertex::Vertex(float x, float y, float z): mHEdge(nullptr), mFlag(0), mKnown(false) {
	mPosition = Eigen::Vector3f(x, y, z);
	mColor = VCOLOR_BLUE;
	mNormal = Eigen::Vector3f::Zero();
}


const Eigen::Vector3f& Vertex::position() const {
	return mPosition;
}

const Eigen::Vector3f& Vertex::setPosition(const Eigen::Vector3f& p) {
	mPosition = p;
	return mPosition;
}

const Eigen::Vector3f& Vertex::normal() const {
	return mNormal;
}

const Eigen::Vector3f& Vertex::setNormal(const Eigen::Vector3f& n) {
	mNormal = n;
	return mNormal;
}

const Eigen::Vector3f& Vertex::color() const {
	return mColor;
}

const Eigen::Vector3f& Vertex::setColor(const Eigen::Vector3f& c) {
	mColor = c;
	return mColor;
}

HEdge* Vertex::halfEdge() const {
	return mHEdge;
}

HEdge* Vertex::setHalfEdge(HEdge* he) {
	mHEdge = he;
	return mHEdge;
}

int Vertex::index() const {
	return mIndex;
}

int Vertex::setIndex(int i) {
	mIndex = i;
	return mIndex;
}

int Vertex::flag() const {
	return mFlag;
}

int Vertex::setFlag(int f) {
	mFlag = f;
	return mFlag;
}

bool Vertex::isValid() const {
	return mValid;
}

bool Vertex::setValid(bool b) {
	mValid = b;
	return mValid;
}

bool Vertex::isKnown() const {
	return mKnown;
}

bool Vertex::setKnown(bool b) {
	mKnown = b;
	return mKnown;
}

bool Vertex::isBoundary() const {
	OneRingHEdge ring(this);
	HEdge* curr = nullptr;
	while (curr = ring.nextHEdge()) {
		if (curr->isBoundary()) {
			return true;
		}
	}
	return false;
}

int Vertex::valence() const {
	int count = 0;
	OneRingVertex ring(this);
	Vertex* curr = nullptr;
	while (curr = ring.nextVertex()) {
		++count;
	}
	return count;
}

Face::Face() : mHEdge(nullptr), mValid(true) {
}

HEdge* Face::halfEdge() const {
	return mHEdge;
}

HEdge* Face::setHalfEdge(HEdge* he) {
	mHEdge = he;
	return mHEdge;
}

bool Face::isBoundary() const {
	HEdge* curr = mHEdge;
	do {
		if (curr->twin()->isBoundary()) {
			return true;
		}
		curr = curr->next();
	} while (curr != mHEdge);
	return false;
}

bool Face::isValid() const {
	return mValid;
}

bool Face::setValid(bool b) {
	mValid = b;
	return mValid;
}

Mesh::Mesh() {
	mVertexPosFlag = true;
	mVertexNormalFlag = true;
	mVertexColorFlag = true;
}

Mesh::~Mesh() {
	clear();
}

const std::vector< HEdge* >& Mesh::edges() const {
	return mHEdgeList;
}

const std::vector< HEdge* >& Mesh::boundaryEdges() const {
	return mBHEdgeList;
}

const std::vector< Vertex* >& Mesh::vertices() const {
	return mVertexList;
}

const std::vector< Face* >& Mesh::faces() const {
	return mFaceList;
}


bool Mesh::isVertexPosDirty() const {
	return mVertexPosFlag;
}

void Mesh::setVertexPosDirty(bool b) {
	mVertexPosFlag = b;
}

bool Mesh::isVertexNormalDirty() const {
	return mVertexNormalFlag;
}

void Mesh::setVertexNormalDirty(bool b) {
	mVertexNormalFlag = b;
}

bool Mesh::isVertexColorDirty() const {
	return mVertexColorFlag;
}

void Mesh::setVertexColorDirty(bool b) {
	mVertexColorFlag = b;
}

bool Mesh::loadMeshFile(const std::string filename) {
	// Use libigl to parse the mesh file
	bool iglFlag = igl::read_triangle_mesh(filename, mVertexMat, mFaceMat);
	if (iglFlag) {
		clear();

		// Construct the half-edge data structure.
		int numVertices = mVertexMat.rows();
		int numFaces = mFaceMat.rows();

		// Fill in the vertex list
		for (int vidx = 0; vidx < numVertices; ++vidx) {
			mVertexList.push_back(new Vertex(mVertexMat(vidx, 0),
			                                 mVertexMat(vidx, 1),
			                                 mVertexMat(vidx, 2)));
		}
		// Fill in the face list
		for (int fidx = 0; fidx < numFaces; ++fidx) {
			addFace(mFaceMat(fidx, 0), mFaceMat(fidx, 1), mFaceMat(fidx, 2));
		}

		std::vector< HEdge* > hedgeList;
		for (int i = 0; i < mBHEdgeList.size(); ++i) {
			if (mBHEdgeList[i]->start()) {
				hedgeList.push_back(mBHEdgeList[i]);
			}
			// TODO
		}
		mBHEdgeList = hedgeList;

		for (int i = 0; i < mVertexList.size(); ++i) {
			mVertexList[i]->adjHEdges.clear();
			mVertexList[i]->setIndex(i);
			mVertexList[i]->setFlag(0);
		}
	} else {
		std::cout << __FUNCTION__ << ": mesh file loading failed!\n";
	}
	return iglFlag;
}

static void _setPrevNext(HEdge* e1, HEdge* e2) {
	e1->setNext(e2);
	e2->setPrev(e1);
}

static void _setTwin(HEdge* e1, HEdge* e2) {
	e1->setTwin(e2);
	e2->setTwin(e1);
}

static void _setFace(Face* f, HEdge* e) {
	f->setHalfEdge(e);
	e->setFace(f);
}

void Mesh::addFace(int v1, int v2, int v3) {
	Face* face = new Face();

	HEdge* hedge[3];
	HEdge* bhedge[3]; // Boundary half-edges
	Vertex* vert[3];

	for (int i = 0; i < 3; ++i) {
		hedge[i] = new HEdge();
		bhedge[i] = new HEdge(true);
	}
	vert[0] = mVertexList[v1];
	vert[1] = mVertexList[v2];
	vert[2] = mVertexList[v3];

	// Connect prev-next pointers
	for (int i = 0; i < 3; ++i) {
		_setPrevNext(hedge[i], hedge[(i + 1) % 3]);
		_setPrevNext(bhedge[i], bhedge[(i + 1) % 3]);
	}

	// Connect twin pointers
	_setTwin(hedge[0], bhedge[0]);
	_setTwin(hedge[1], bhedge[2]);
	_setTwin(hedge[2], bhedge[1]);

	// Connect start pointers for bhedge
	bhedge[0]->setStart(vert[1]);
	bhedge[1]->setStart(vert[0]);
	bhedge[2]->setStart(vert[2]);
	for (int i = 0; i < 3; ++i) {
		hedge[i]->setStart(vert[i]);
	}

	// Connect start pointers
	// Connect face-hedge pointers
	for (int i = 0; i < 3; ++i) {
		vert[i]->setHalfEdge(hedge[i]);
		vert[i]->adjHEdges.push_back(hedge[i]);
		_setFace(face, hedge[i]);
	}
	vert[0]->adjHEdges.push_back(bhedge[1]);
	vert[1]->adjHEdges.push_back(bhedge[0]);
	vert[2]->adjHEdges.push_back(bhedge[2]);

	// Merge boundary if needed
	for (int i = 0; i < 3; ++i) {
		Vertex* start = bhedge[i]->start();
		Vertex* end = bhedge[i]->end();

		for (int j = 0; j < end->adjHEdges.size(); ++j) {
			HEdge* curr = end->adjHEdges[j];
			if (curr->isBoundary() && curr->end() == start) {
				_setPrevNext(bhedge[i]->prev(), curr->next());
				_setPrevNext(curr->prev(), bhedge[i]->next());
				_setTwin(bhedge[i]->twin(), curr->twin());
				bhedge[i]->setStart(nullptr); // Mark as unused
				curr->setStart(nullptr); // Mark as unused
				break;
			}
		}
	}

	// Finally add hedges and faces to list
	for (int i = 0; i < 3; ++i) {
		mHEdgeList.push_back(hedge[i]);
		mBHEdgeList.push_back(bhedge[i]);
	}
	mFaceList.push_back(face);
}

Eigen::Vector3f Mesh::initBboxMin() const {
	return (mVertexMat.colwise().minCoeff()).transpose();
}

Eigen::Vector3f Mesh::initBboxMax() const {
	return (mVertexMat.colwise().maxCoeff()).transpose();
}

void Mesh::groupingVertexFlags() {
	// Init to 255
	for (Vertex* vert : mVertexList) {
		if (vert->flag() != 0) {
			vert->setFlag(255);
		}
	}
	// Group handles
	int id = 0;
	std::vector< Vertex* > tmpList;
	for (Vertex* vert : mVertexList) {
		if (vert->flag() == 255) {
			++id;
			vert->setFlag(id);

			// Do search
			tmpList.push_back(vert);
			while (!tmpList.empty()) {
				Vertex* v = tmpList.back();
				tmpList.pop_back();

				OneRingVertex orv = OneRingVertex(v);
				while (Vertex* v2 = orv.nextVertex()) {
					if (v2->flag() == 255) {
						v2->setFlag(id);
						tmpList.push_back(v2);
					}
				}
			}
		}
	}
}

void Mesh::clear() {
	for (int i = 0; i < mHEdgeList.size(); ++i) {
		delete mHEdgeList[i];
	}
	for (int i = 0; i < mBHEdgeList.size(); ++i) {
		delete mBHEdgeList[i];
	}
	for (int i = 0; i < mVertexList.size(); ++i) {
		delete mVertexList[i];
	}
	for (int i = 0; i < mFaceList.size(); ++i) {
		delete mFaceList[i];
	}

	mHEdgeList.clear();
	mBHEdgeList.clear();
	mVertexList.clear();
	mFaceList.clear();
}

std::vector< int > Mesh::collectMeshStats() {
	int V = 0; // # of vertices
	int E = 0; // # of half-edges
	int F = 0; // # of faces
	int B = 0; // # of boundary loops
	int C = 0; // # of connected components
	int G = 0; // # of genus

	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Collect mesh information as listed above.
	/**********************************************/
	V = mVertexList.size();
	F = mFaceList.size();
	E = mHEdgeList.size() / 2 + mBHEdgeList.size();
	C = countConnectedComponents();
	B = countBoundaryLoops();
	G = (2 * C - B - F + E - V) / 2;

	/*====== Programming Assignment 0 ======*/

	std::vector< int > stats;
	stats.push_back(V);
	stats.push_back(E);
	stats.push_back(F);
	stats.push_back(B);
	stats.push_back(C);
	stats.push_back(G);
	return stats;
}


int Mesh::countBoundaryLoops() {
	int count = 0;

	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Helper function for Mesh::collectMeshStats()
	/**********************************************/
	for (int i = 0; i < mBHEdgeList.size(); i++) {
		HEdge* curr = mBHEdgeList[i];
		if (!curr->flag()) {
			curr->setFlag(true);
			Vertex* start = curr->start();
			Vertex* end = curr->end();
			
			int test = 0;
			while(curr->next() != mBHEdgeList[i]){
				curr = curr->next();
				curr->setFlag(true);
            }
			count++;
		}
	}


	/*====== Programming Assignment 0 ======*/

	return count;
}

void DFS(Vertex* v) {
	Vertex* w;
	v->setKnown(true);
	OneRingVertex ring(v);
	while (w = ring.nextVertex()) {
		if (!w->isKnown()) {
			DFS(w);
		}
	}
}

int Mesh::countConnectedComponents() {
	int count = 0;

	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Helper function for Mesh::collectMeshStats()
	/* Count the number of connected components of
	/* the mesh. (Hint: use a stack)
	/**********************************************/
	for (int i = 0; i < mVertexList.size(); i++) {
		if (!mVertexList[i]->isKnown()) {
			DFS(mVertexList[i]);
			count++;
		}
	}
	/*====== Programming Assignment 0 ======*/

	return count;
}


void Mesh::computeVertexNormals() {
	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Compute per-vertex normal using neighboring
	/* facet information. (Hint: remember using a
	/* weighting scheme. Plus, do you notice any
	/* disadvantages of your weighting scheme?)
	/**********************************************/
	for (int i = 0; i < mVertexList.size(); i++) {
		Vertex* v = mVertexList[i];
		int cnt = v->valence();
		Eigen::Vector3f posV = v->position();
		Eigen::Vector3f sum = Eigen::Vector3f::Zero();
		OneRingHEdge ring(v);
		HEdge* curr = ring.nextHEdge();
		Eigen::Vector3f posCurr = curr->end()->position();
		Eigen::Vector3f vecCurr = posCurr - posV;
		while (HEdge* next = ring.nextHEdge()) {
			Eigen::Vector3f posNext = next->end()->position();
			Eigen::Vector3f vecNext = posNext - posV;
			Eigen::Vector3f multi = vecCurr.cross(vecNext);
			sum += multi / cnt;
			vecCurr = vecNext;
		}
		v->setNormal(sum);
	}


	/*====== Programming Assignment 0 ======*/

	// Notify mesh shaders
	setVertexNormalDirty(true);
}


void Mesh::umbrellaSmooth(bool cotangentWeights) {
	/*====== Programming Assignment 1 ======*/
	float lamda = 0.9;
	//construst the matrix of all the vectors
	Eigen::MatrixXd matX = Eigen::MatrixXd::Zero(mVertexList.size(), 3);
	for (int i = 0; i < mVertexList.size(); i++) {
		Eigen::Vector3d vec = mVertexList[i]->position().cast<double>();
		matX.row(i) = vec.transpose();
	}
	//std::cout << mVertexList[0]->position()<<std::endl<<std::endl;
	//std::cout << matX.row(0) << std::endl << std::endl;

	if (cotangentWeights) {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 1: Implement the cotangent weighting
		/* scheme for explicit mesh smoothing.
		/*
		/* Hint:
		/* It is advised to double type to store the
		/* weights to avoid numerical issues.
		/**********************************************/
		//initial L matrix
		Eigen::SparseMatrix<double> mat(mVertexList.size(), mVertexList.size());
		for (int i = 0; i < mVertexList.size(); i++)
			mat.coeffRef(i, i) = -1;
		//std::cout << mat.nonZeros() << std::endl << std::endl;
		//compute the L matrix
		for (int i = 0; i < mVertexList.size(); i++) {
			Vertex* core = mVertexList[i];
			std::vector< Vertex* > adjVertexList;
			OneRingVertex ring(mVertexList[i]);
			while (Vertex* v = ring.nextVertex()) {
				adjVertexList.push_back(v);
			}
			double sumW = 0;
			for (int j = 0; j < core->valence(); j++) {
				Vertex* curr = adjVertexList[j], *next = nullptr, *prev = nullptr;
				if (j == 0) {
					prev = adjVertexList[core->valence() - 1];
					next = adjVertexList[j + 1];
				}
				else if (j == core->valence() - 1) {
					next = adjVertexList[0];
					prev = adjVertexList[j - 1];
				}
				else {
					prev = adjVertexList[j - 1];
					next = adjVertexList[j + 1];
				}
				Eigen::Vector3f vecAlpha1 = core->position() - next->position();
				Eigen::Vector3f vecAlpha2 = curr->position() - next->position();

				double cosAlpha = vecAlpha1.dot(vecAlpha2) / sqrt(vecAlpha1.dot(vecAlpha1) * vecAlpha2.dot(vecAlpha2));
				double cotAlpha = cosAlpha / sin(acos(cosAlpha));

				Eigen::Vector3f vecBeta1 = core->position() - prev->position();
				Eigen::Vector3f vecBeta2 = curr->position() - prev->position();
				double cosBeta = vecBeta1.dot(vecBeta2) / sqrt(vecBeta1.dot(vecBeta1) * vecBeta2.dot(vecBeta2));
				double cotBeta = cosBeta / sin(acos(cosBeta));

				double w = cotAlpha + cotBeta;
				mat.coeffRef(core->index(), curr->index()) = w;
				sumW += w;
			}
			for (int j = 0; j < core->valence(); j++){
				Vertex* curr = adjVertexList[j];
				mat.coeffRef(core->index(), curr->index()) = mat.coeffRef(core->index(), curr->index()) / sumW;
			}
		}
		//Xt+1 = Xt + lamda * L * Xt
		Eigen::MatrixXd matRes = Eigen::MatrixXd::Zero(mVertexList.size(), 3);
		matRes = matX + lamda * mat * matX;
		for (int i = 0; i < mVertexList.size(); i++) {
			Vertex* core = mVertexList[i];
			Eigen::Vector3d updatePos = matRes.row(i).transpose();
			core->setPosition(updatePos.cast<float>());
		}



	} else {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 2: Implement the uniform weighting
		/* scheme for explicit mesh smoothing.
		/**********************************************/
		for (int i = 0; i < mVertexList.size(); i++) {
			Vertex* core = mVertexList[i];
			Eigen::Vector3f uniformLP = Eigen::Vector3f::Zero();
			OneRingVertex ring(core);
			Vertex* adj = nullptr;
			while (adj = ring.nextVertex()) {
				uniformLP += adj->position();
			}
			uniformLP = uniformLP / core->valence() - core->position();
			core->setPosition(core->position() + lamda * uniformLP);
		}
	}

	/*====== Programming Assignment 1 ======*/

	computeVertexNormals();
	// Notify mesh shaders
	setVertexPosDirty(true);
}

void Mesh::implicitUmbrellaSmooth(bool cotangentWeights) {
	/*====== Programming Assignment 1 ======*/
	float lamda = 1.0;

	Eigen::VectorXd bx(mVertexList.size()), by(mVertexList.size()), bz(mVertexList.size()), x_sol(mVertexList.size()), y_sol(mVertexList.size()), z_sol(mVertexList.size());
	for (int k = 0; k < mVertexList.size(); k++) {
		bx[k] = mVertexList[k]->position()[0];
		by[k] = mVertexList[k]->position()[1];
		bz[k] = mVertexList[k]->position()[2];
		x_sol[k] = y_sol[k] = z_sol[k] = 0;
	}



	/* A sparse linear system Ax=b solver using the conjugate gradient method. */
	auto fnConjugateGradient = [](const Eigen::SparseMatrix< double >& A,
	                              const Eigen::VectorXd& b,
	                              int maxIterations,
	                              float errorTolerance,
	                              Eigen::VectorXd& x)
	{
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Params:
		/*  A:
		/*  b:
		/*  maxIterations:	Max number of iterations
		/*  errorTolerance: Error tolerance for the early stopping condition
		/*  x:				Stores the final solution, but should be initialized.
		/**********************************************/
		/*
		/* Step 1: Implement the conjugate gradient
		/* method.
		/* Hint: https://en.wikipedia.org/wiki/Conjugate_gradient_method
		/**********************************************/
		/*
		// conjugater gradient
		int i = 0;
		double error = 100;
		Eigen::VectorXd x0(A.rows());
		x0.fill(1);
		Eigen::VectorXd r0(A.rows());
		Eigen::VectorXd r0c(A.rows());

		Eigen::VectorXd r1(A.rows());
		Eigen::VectorXd r1c(A.rows());

		Eigen::VectorXd p0(A.rows());
		Eigen::VectorXd p0c(A.rows());

		Eigen::VectorXd p1(A.rows());
		Eigen::VectorXd p1c(A.rows());

		r0 = b - A*x0;
		//r0c=b.conjugate()-x0.conjugate()*A.transpose();
		p0 = r0;

		while (i<maxIterations&&error>errorTolerance) {
			double a0;
			double b0;

			a0 = double(r0.transpose()*r0) / double(p0.transpose()*A*p0);
			x = x0 + a0*p0;
			r1 = r0 - a0*A*p0;

			b0 = double(r1.transpose()*r1) / double(r0.transpose()*r0);
			p1 = r1 + b0*p0;

			p0 = p1;
			r0 = r1;
			x0 = x;
			i++;
			error = r1.squaredNorm();
		}*/

		//biconjudgate gradient 
		Eigen::VectorXd r, rr, r_next, p, v, p_next, v_next, h, x_next, s, t;
		r = b - A * x;
		rr = r;
		double pp = 1, alpha = 1, w = 1, pp_next, beta, w_next;
		p = Eigen::VectorXd::Zero(A.rows());
		v = Eigen::VectorXd::Zero(A.rows());
		for (int i = 1; i < maxIterations; i++) {
			pp_next = rr.transpose() * r;
			beta = (pp_next / pp) * (alpha / w);
			p_next = r + beta * (p - w * v);
			v_next = A * p_next;
			alpha = pp_next / (rr.transpose() * v_next);
			h = x + alpha * p_next;
			if (h.transpose() * h < errorTolerance) {
				x = h;
				break;
			}
			s = r - alpha * v_next;
			t = A * s;
			double tmp1 = t.transpose() * s;
			double tmp2 = t.transpose() * t;
			w_next = tmp1 / tmp2;
			x_next = h + w_next * s;
			if (x_next.transpose() * x_next < errorTolerance) {
				x = x_next;
				break;
			}
			r_next = s - w_next * t;
			r = r_next; p = p_next; v = v_next; x = x_next; pp = pp_next; w = w_next;
		}



	};

	/* IMPORTANT:
	/* Please refer to the following link about the sparse matrix construction in Eigen. */
	/* http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title3 */

	//initial L matrix
	Eigen::SparseMatrix<double> mat(mVertexList.size(), mVertexList.size());
	for (int i = 0; i < mVertexList.size(); i++)
		mat.coeffRef(i, i) = -1;

	//initial an identify matrix
	Eigen::SparseMatrix<double> id(mVertexList.size(), mVertexList.size());
	for (int i = 0; i < mVertexList.size(); i++)
		id.coeffRef(i, i) = 1;

	Eigen::SparseMatrix<double> matRes(mVertexList.size(), mVertexList.size());
	if (cotangentWeights) {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 2: Implement the cotangent weighting
		/* scheme for implicit mesh smoothing. Use
		/* the above fnConjugateGradient for solving
		/* sparse linear systems.
		/*
		/* Hint:
		/* It is advised to double type to store the
		/* weights to avoid numerical issues.
		/**********************************************/
		//compute the L matrix
		for (int i = 0; i < mVertexList.size(); i++) {
			Vertex* core = mVertexList[i];
			std::vector< Vertex* > adjVertexList;
			OneRingVertex ring(mVertexList[i]);
			while (Vertex* v = ring.nextVertex()) {
				adjVertexList.push_back(v);
			}
			double sumW = 0;
			for (int j = 0; j < core->valence(); j++) {
				Vertex* curr = adjVertexList[j], *next = nullptr, *prev = nullptr;
				if (j == 0) {
					prev = adjVertexList[core->valence() - 1];
					next = adjVertexList[j + 1];
				}
				else if (j == core->valence() - 1) {
					next = adjVertexList[0];
					prev = adjVertexList[j - 1];
				}
				else {
					prev = adjVertexList[j - 1];
					next = adjVertexList[j + 1];
				}
				Eigen::Vector3f vecAlpha1 = core->position() - next->position();
				Eigen::Vector3f vecAlpha2 = curr->position() - next->position();

				double cosAlpha = vecAlpha1.dot(vecAlpha2) / sqrt(vecAlpha1.dot(vecAlpha1) * vecAlpha2.dot(vecAlpha2));
				double cotAlpha = cosAlpha / sin(acos(cosAlpha));

				Eigen::Vector3f vecBeta1 = core->position() - prev->position();
				Eigen::Vector3f vecBeta2 = curr->position() - prev->position();
				double cosBeta = vecBeta1.dot(vecBeta2) / sqrt(vecBeta1.dot(vecBeta1) * vecBeta2.dot(vecBeta2));
				double cotBeta = cosBeta / sin(acos(cosBeta));

				double w = cotAlpha + cotBeta;
				mat.coeffRef(core->index(), curr->index()) = w;
				sumW += w;
			}
			for (int j = 0; j < core->valence(); j++) {
				Vertex* curr = adjVertexList[j];
				mat.coeffRef(core->index(), curr->index()) = mat.coeffRef(core->index(), curr->index()) / sumW;
			}
		}
		// E - lamda * L
		
		matRes = id - lamda * mat;
		


	} else {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 3: Implement the uniform weighting
		/* scheme for implicit mesh smoothing. Use
		/* the above fnConjugateGradient for solving
		/* sparse linear systems.
		/**********************************************/

		//compute L matrix
		for (int i = 0; i < mVertexList.size(); i++) {
			Vertex* core = mVertexList[i];
			std::vector< Vertex* > adjVertexList;
			OneRingVertex ring(mVertexList[i]);
			while (Vertex* v = ring.nextVertex()) {
				adjVertexList.push_back(v);
			}
			double sumW = 0;
			for (int j = 0; j < core->valence(); j++) {
				Vertex* curr = adjVertexList[j];				
				double w = 1;
				mat.coeffRef(core->index(), curr->index()) = w;
				sumW += w;
			}
			for (int j = 0; j < core->valence(); j++) {
				Vertex* curr = adjVertexList[j];
				mat.coeffRef(core->index(), curr->index()) = mat.coeffRef(core->index(), curr->index()) / sumW;
			}
		}
		// E - lamda * L

		matRes = id - lamda * mat;

	}

	//call the conjugate gradient
	fnConjugateGradient(matRes, bx, 100, 1e-6, x_sol);
	fnConjugateGradient(matRes, by, 100, 1e-6, y_sol);
	fnConjugateGradient(matRes, bz, 100, 1e-6, z_sol);
	for (int k = 0; k < mVertexList.size(); k++) {
		//			std::cout << x_sol[k] << std::endl;
		mVertexList[k]->setPosition(Eigen::Vector3f((float)x_sol[k], (float)y_sol[k], (float)z_sol[k]));
	}
	/*====== Programming Assignment 1 ======*/

	computeVertexNormals();
	// Notify mesh shaders
	setVertexPosDirty(true);
}
