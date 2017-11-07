#ifndef DEFORMER_H
#define DEFORMER_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mesh.h"

// Deform mesh using Laplacian coordinates
class Deformer {
public:
	Deformer();
	~Deformer();

	void setMesh(Mesh* mesh);
	//void setRotInva(bool mRotInva);

	/*====== Programming Assignment 2 ======*/
	// This is the place where the editing techniques take place
	void deform();
	/*====== Programming Assignment 2 ======*/

private:
	/*====== Programming Assignment 2 ======*/
	// Build left hand side matrix and pre-factorize it
	void buildSystemMat();
	/*====== Programming Assignment 2 ======*/

	void clear();

	Mesh* mMesh;
	std::vector< Vertex* > mRoiList;
	// Solver for pre-factorizing the system matrix of the deformation
	Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >* mCholeskySolver;
	// sparse matrices of A, A.transpose(), and B.
	Eigen::SparseMatrix<double> smA;
	Eigen::SparseMatrix<double> smAT;
	Eigen::SparseMatrix<double> smB;
	//flag for Laplacian editing whether it's ratiation invariant or not
	//initial as false for naive Laplacian editing
	bool rotation_invariant;
};

#endif // DEFORMER_H
