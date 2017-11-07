#include "deformer.h"
#include <iostream>

Deformer::Deformer() : mMesh(nullptr),
				       rotation_invariant(false),  //initial as false for naive Laplacian editing
                       mCholeskySolver(nullptr) {
}

Deformer::~Deformer() {
	clear();
}

void Deformer::clear() {
	if (mCholeskySolver) {
		delete mCholeskySolver;
	}
	mCholeskySolver = nullptr;
	mRoiList.clear();
}


void Deformer::setMesh(Mesh* mesh) {
	mMesh = mesh;
	clear();
	// Record the handle vertices
	for (Vertex* vert : mMesh->vertices()) {
		if (vert->flag() > 0 || vert->isBoundary()) {
			mRoiList.push_back(vert);
		}
	}
	// Build system matrix for deformation
	buildSystemMat();
}


//construct sparse matrices
Eigen::SparseMatrix<double> toSparseMatrix(int rows, int cols, std::vector< Eigen::Triplet<double>> tmpV) {
	Eigen::SparseMatrix<double> tmp(rows, cols);
	tmp.setFromTriplets(tmpV.begin(), tmpV.end());
	return tmp;
}

void Deformer::buildSystemMat() {
	/*====== Programming Assignment 2 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Build the matrix of the linear system for 
	/* deformation and do factorization, in order
	/* to reuse and speed up in Deformer::deform().
	/* Handle vertices are maked by Vertex::flag() > 0
	/* Movements of the specified handle are already
	/* recorded in Vertex::position()
	/**********************************************/

	Eigen::SparseMatrix< double > systemMat;

	/*====== Programming Assignment 2 ======*/

	// Please refer to the following link for the usage of sparse linear system solvers in Eigen
	// https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

	//used for construct the sparse matrix
	// for naive laplacian editing
	std::vector< Eigen::Triplet<double>> L;
	std::vector< Eigen::Triplet<double>> A; 
	std::vector< Eigen::Triplet<double>> B;
	// for laplacian editing with local rotations
	std::vector< Eigen::Triplet<double>> T;
	std::vector< Eigen::Triplet<double>> A_r; 
	std::vector< Eigen::Triplet<double>> B_r;
	
	std::vector< Vertex* > vList = mMesh->vertices();
	int N = vList.size();

	// 1. compute cotangent weights and construct the Laplacian matrix
	for (Vertex* vert : mMesh->vertices()) {
		OneRingVertex ring(vert); 
		std::vector< Vertex* > adjVertexList;
		std::vector< Eigen::Triplet<double>> vWeight(0);
		
		double sumW = 0;
		while (Vertex* v = ring.nextVertex()) {
			adjVertexList.push_back(v);
		}
		for (int j = 0; j < vert->valence(); j++) {
			Vertex* curr = adjVertexList[j], *next = nullptr, *prev = nullptr;
			if (j == 0) {
				prev = adjVertexList[vert->valence() - 1];
				next = adjVertexList[j + 1];
			}
			else if (j == vert->valence() - 1) {
				next = adjVertexList[0];
				prev = adjVertexList[j - 1];
			}
			else {
				prev = adjVertexList[j - 1];
				next = adjVertexList[j + 1];
			}
			Eigen::Vector3f vecAlpha1 = vert->position() - next->position();
			Eigen::Vector3f vecAlpha2 = curr->position() - next->position();

			double cosAlpha = vecAlpha1.dot(vecAlpha2) / sqrt(vecAlpha1.dot(vecAlpha1) * vecAlpha2.dot(vecAlpha2));
			double cotAlpha = cosAlpha / sin(acos(cosAlpha));

			Eigen::Vector3f vecBeta1 = vert->position() - prev->position();
			Eigen::Vector3f vecBeta2 = curr->position() - prev->position();
			double cosBeta = vecBeta1.dot(vecBeta2) / sqrt(vecBeta1.dot(vecBeta1) * vecBeta2.dot(vecBeta2));
			double cotBeta = cosBeta / sin(acos(cosBeta));

			double w = cotAlpha + cotBeta;
			sumW += w;
			
			vWeight.push_back(Eigen::Triplet<double>(vert->index(), curr->index(), w)); 
		}
		for (int j = 0; j < vWeight.size(); j++) {
			//(i,j,w) -- (i, j) is the position of the weight in the Laplacian matrix
			L.push_back(Eigen::Triplet<double>(vert->index(), vWeight[j].col(), -vWeight[j].value() / sumW));
			A.push_back(Eigen::Triplet<double>(vert->index(), vWeight[j].col(), -vWeight[j].value() / sumW));
			A_r.push_back(Eigen::Triplet<double>(vert->index(), vWeight[j].col(), -vWeight[j].value() / sumW));//x
			A_r.push_back(Eigen::Triplet<double>(vert->index() + N, vWeight[j].col() + N, -vWeight[j].value() / sumW));//y
			A_r.push_back(Eigen::Triplet<double>(vert->index() + N * 2, vWeight[j].col() + N * 2, -vWeight[j].value() / sumW));//z
		}
		//(i,i,1)
		L.push_back(Eigen::Triplet<double>(vert->index(), vert->index(), 1));
		A.push_back(Eigen::Triplet<double>(vert->index(), vert->index(), 1));
		A_r.push_back(Eigen::Triplet<double>(vert->index(), vert->index(), 1)); //x
		A_r.push_back(Eigen::Triplet<double>(vert->index() + N, vert->index() + N, 1)); //y
		A_r.push_back(Eigen::Triplet<double>(vert->index() + 2 * N, vert->index() + 2 * N, 1));//z
	}
	
	Eigen::SparseMatrix<double> smL(vList.size(), vList.size());
	smL.setFromTriplets(L.begin(), L.end());    // Lapcian matrix


	//2. compute original Laplacian coordinates
	// construct vertex matrix
	std::vector< Eigen::Triplet<double>> V;
	for (Vertex* vert : mMesh->vertices()) {
		Eigen::Vector3d pos = vert->position().cast<double>();
		V.push_back(Eigen::Triplet<double>(vert->index(), 0, pos[0]));
		V.push_back(Eigen::Triplet<double>(vert->index(), 1, pos[1]));
		V.push_back(Eigen::Triplet<double>(vert->index(), 2, pos[2]));
	}
	Eigen::SparseMatrix<double> smV(vList.size(), 3);
	smV.setFromTriplets(V.begin(), V.end());
	//store the original Laplacian coordinates in orgLap
	Eigen::SparseMatrix<double> orgLap = smL * smV;


	// 3. build systemMat
	for (int i = 0; i < mRoiList.size(); i++) {
		Vertex* conVertex = mRoiList[i];
		//Eigen::Vector3d pos = conVertex->position().cast<double>();
		int index = conVertex->index();
		A.push_back(Eigen::Triplet<double>(i + vList.size(), index, 1));
		A_r.push_back(Eigen::Triplet<double>(i + 3 * N, index, 1)); 
		A_r.push_back(Eigen::Triplet<double>(i + mRoiList.size() + 3 * N, index + N, 1));
		A_r.push_back(Eigen::Triplet<double>(i + mRoiList.size() * 2 + 3 * N, index + 2 * N, 1));
	}
	
	

	// naive Laplacian editing
	if (rotation_invariant == false) {
		//construct the sparse matrix A 
		smA = toSparseMatrix(vList.size() + mRoiList.size(), vList.size(), A);

		for (int i = 0; i < vList.size(); i++) {
			B.push_back(Eigen::Triplet<double>(i, 0, orgLap.coeff(i, 0)));
			B.push_back(Eigen::Triplet<double>(i, 1, orgLap.coeff(i, 1)));
			B.push_back(Eigen::Triplet<double>(i, 2, orgLap.coeff(i, 2)));
		}
		//compute B
		//Eigen::SparseMatrix<double> smB(vList.size() + mRoiList.size(), 3);
		//smB.setFromTriplets(B.begin(), B.end());
		smB = toSparseMatrix(vList.size() + mRoiList.size(), 3, B);
	}
	// Lapacian edting with local transition
	// the size of the solution V' is 3N * 1 (x, y, z on the same column) 
	else {
		//construct the sparse matrix A
		smA = toSparseMatrix(3 * (vList.size() + mRoiList.size()), 3 * vList.size(), A_r);

		for (int i = 0; i < vList.size(); i++) {
			Vertex* vert = vList[i];
			int cnt = vert->valence();
			// construct Ai matrix for T
			Eigen::MatrixXd Ai(3 * (cnt + 1), 7);
			Ai(0, 0) = vert->position()[0];
			Ai(0, 1) = 0;
			Ai(0, 2) = vert->position()[2];
			Ai(0, 3) = -vert->position()[1];
			Ai(0, 4) = 1;
			Ai(0, 5) = 0;
			Ai(0, 6) = 0;

			Ai(1, 0) = vert->position()[1];
			Ai(1, 1) = -vert->position()[2];
			Ai(1, 2) = 0;
			Ai(1, 3) = vert->position()[0];
			Ai(1, 4) = 0;
			Ai(1, 5) = 1;
			Ai(1, 6) = 0;

			Ai(2, 0) = vert->position()[2];
			Ai(2, 1) = vert->position()[1];
			Ai(2, 2) = -vert->position()[0];
			Ai(2, 3) = 0;
			Ai(2, 4) = 0;
			Ai(2, 5) = 0;
			Ai(2, 6) = 1;

			OneRingVertex ring(vert);
			for (int j = 1; j < cnt+1; j++) {
				Vertex* neighbor = ring.nextVertex();
				Ai(j * 3, 0) = neighbor->position()[0];
				Ai(j * 3, 1) = 0;
				Ai(j * 3, 2) = neighbor->position()[2];
				Ai(j * 3, 3) = -neighbor->position()[1];
				Ai(j * 3, 4) = 1;
				Ai(j * 3, 5) = 0;
				Ai(j * 3, 6) = 0;

				Ai(j * 3 + 1, 0) = neighbor->position()[1];
				Ai(j * 3 + 1, 1) = -neighbor->position()[2];
				Ai(j * 3 + 1, 2) = 0;
				Ai(j * 3 + 1, 3) = neighbor->position()[0];
				Ai(j * 3 + 1, 4) = 0;
				Ai(j * 3 + 1, 5) = 1;
				Ai(j * 3 + 1, 6) = 0;

				Ai(j * 3 + 2, 0) = neighbor->position()[2];
				Ai(j * 3 + 2, 1) = neighbor->position()[1];
				Ai(j * 3 + 2, 2) = -neighbor->position()[0];
				Ai(j * 3 + 2, 3) = 0;
				Ai(j * 3 + 2, 4) = 0;
				Ai(j * 3 + 2, 5) = 0;
				Ai(j * 3 + 2, 6) = 1;
			}
			// compute the pseudo inverse of Ai
			Eigen::MatrixXd invAi = eigenPinv(Ai);
			// construct Ti * original Laplacian coordinates
			double xLap = orgLap.coeff(i, 0);
			double yLap = orgLap.coeff(i, 1);
			double zLap = orgLap.coeff(i, 2);
			//editted Laplacian coordinate x
			// s * x - h3 * y + h2 * z
			T.push_back(Eigen::Triplet<double>(i, i, xLap * invAi(0, 0) - yLap * invAi(3, 0) + zLap * invAi(2, 0))); 
			T.push_back(Eigen::Triplet<double>(i, i+N, xLap * invAi(0, 1) - yLap * invAi(3, 1) + zLap * invAi(2, 1)));
			T.push_back(Eigen::Triplet<double>(i, i+2*N, xLap * invAi(0, 2) - yLap * invAi(3, 2) + zLap * invAi(2, 2))); 
			//editted Laplacian coordinate y
			//h3 * x + s * y - h1 * z
			T.push_back(Eigen::Triplet<double>(i + N, i, xLap * invAi(3, 0) + yLap * invAi(0, 0) - zLap * invAi(1, 0)));
			T.push_back(Eigen::Triplet<double>(i + N, i+N, xLap * invAi(3, 1) + yLap * invAi(0, 1) - zLap * invAi(1, 1)));
			T.push_back(Eigen::Triplet<double>(i + N, i+2*N, xLap * invAi(3, 2) + yLap * invAi(0, 2) - zLap * invAi(1, 2)));
			//editted Laplacian coordinate z
			//-h2 * x + h1 * y+s * z
			T.push_back(Eigen::Triplet<double>(i + N * 2, i, -xLap * invAi(2, 0) + yLap * invAi(1, 0) + zLap * invAi(0, 0)));
			T.push_back(Eigen::Triplet<double>(i + N * 2, i+N, -xLap * invAi(2, 1) + yLap * invAi(1, 1) + zLap * invAi(0, 1)));
			T.push_back(Eigen::Triplet<double>(i + N * 2, i+2*N, -xLap * invAi(2, 2) + yLap * invAi(1, 2) + zLap * invAi(0, 2)));

			OneRingVertex ring2(vert);
			for (int j = 1; j < cnt+1; j++) {
				Vertex* neighbor = ring2.nextVertex();
				int index = neighbor->index();
				T.push_back(Eigen::Triplet<double>(i, index, xLap * invAi(0, 3*j) - yLap * invAi(3, 3*j) + zLap * invAi(2, 3*j)));
				T.push_back(Eigen::Triplet<double>(i, index + N, xLap * invAi(0, 3*j+1) - yLap * invAi(3, 3 * j + 1) + zLap * invAi(2, 3 * j + 1)));
				T.push_back(Eigen::Triplet<double>(i, index + 2 * N, xLap * invAi(0, 3 * j + 2) - yLap * invAi(3, 3 * j + 2) + zLap * invAi(2, 3 * j + 2)));

				T.push_back(Eigen::Triplet<double>(i + N, index, xLap * invAi(3, 3 * j) + yLap * invAi(0, 3 * j) - zLap * invAi(1, 3 * j)));
				T.push_back(Eigen::Triplet<double>(i + N, index + N, xLap * invAi(3, 3 * j + 1) + yLap * invAi(0, 3 * j + 1) - zLap * invAi(1, 3 * j + 1)));
				T.push_back(Eigen::Triplet<double>(i + N, index + 2 * N, xLap * invAi(3, 3 * j + 2) + yLap * invAi(0, 3 * j + 2) - zLap * invAi(1, 3 * j + 2)));

				T.push_back(Eigen::Triplet<double>(i + N * 2, index, -xLap * invAi(2, 3 * j) + yLap * invAi(1, 3 * j) + zLap * invAi(0, 3 * j)));
				T.push_back(Eigen::Triplet<double>(i + N * 2, index + N, -xLap * invAi(2, 3 * j + 1) + yLap * invAi(1, 3 * j + 1) + zLap * invAi(0, 3 * j + 1)));
				T.push_back(Eigen::Triplet<double>(i + N * 2, index + 2 * N, -xLap * invAi(2, 3 * j + 2) + yLap * invAi(1, 3 * j + 2) + zLap * invAi(0, 3 * j + 2)));

			}
		} // end-i
		Eigen::SparseMatrix<double> smT = toSparseMatrix(3 * (N + mRoiList.size()), 3 * N, T);
		smA = smA - smT;

		smB = toSparseMatrix(3 * (N + mRoiList.size()), 1, B);
	}

	// compute the system matrix
	smAT = smA.transpose();
	systemMat = smAT * smA;
	

	// Do factorization
	if (systemMat.nonZeros() > 0) {
		mCholeskySolver = new Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >();
		mCholeskySolver->compute(systemMat);
		if (mCholeskySolver->info() != Eigen::Success) {
			// Decomposition failed
			std::cout << "Sparse decomposition failed\n";
		} else {
			std::cout << "Sparse decomposition succeeded\n";
		} 
	}
}



void Deformer::deform() {
	if (mCholeskySolver == nullptr) {
		return;
	}

	/*====== Programming Assignment 2 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* This is the place where the editing techniques 
	/* take place.
	/* Solve for the new vertex positions after the 
	/* specified handles move using the factorized
	/* matrix from Deformer::buildSystemMat(), i.e.,
	/* mCholeskySolver defined in deformer.h
	/**********************************************/

	// Please refer to the following link for the usage of sparse linear system solvers in Eigen
	// https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

	std::vector< Vertex* > vList = mMesh->vertices();
	//naive Laplacian editing
	if (rotation_invariant == false) {
		for (int i = 0; i < mRoiList.size(); i++) {
			Eigen::Vector3d pos = mRoiList[i]->position().cast<double>();
			smB.coeffRef(i + vList.size(), 0) = pos[0];
			smB.coeffRef(i + vList.size(), 1) = pos[1];
			smB.coeffRef(i + vList.size(), 2) = pos[2];
		}
		Eigen::SparseMatrix<double> smAT_B = smAT * smB;

		Eigen::SparseMatrix<double> newV = mCholeskySolver->solve(smAT_B);
		for (Vertex* vert : mMesh->vertices()) {
			Eigen::Vector3f v;
			v[0] = (float)newV.coeff(vert->index(), 0);
			v[1] = (float)newV.coeff(vert->index(), 1);
			v[2] = (float)newV.coeff(vert->index(), 2);
			vert->setPosition(v);
		}
	}
	//Laplacian editing with local rotations
	else {
		for (int i = 0; i < mRoiList.size(); i++) {
			Eigen::Vector3d pos = mRoiList[i]->position().cast<double>();
			smB.coeffRef(i + vList.size() * 3, 0) = pos[0];
			smB.coeffRef(i + mRoiList.size() + vList.size() * 3, 0) = pos[1];
			smB.coeffRef(i + mRoiList.size() * 2 + vList.size() * 3, 0) = pos[2];
		}
		Eigen::SparseMatrix<double> smAT_B = smAT * smB;

		Eigen::SparseMatrix<double> newV = mCholeskySolver->solve(smAT_B);
		for (Vertex* vert : mMesh->vertices()) {
			Eigen::Vector3f v;
			v[0] = (float)newV.coeff(vert->index(), 0);
			v[1] = (float)newV.coeff(vert->index() + vList.size(), 0);
			v[2] = (float)newV.coeff(vert->index() + vList.size() * 2, 0);
			vert->setPosition(v);
		}
	}

	/*====== Programming Assignment 2 ======*/
}
