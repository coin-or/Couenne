/* $Id$
 *
 * Name:    orthocut.cpp
 * Author:  Andrea Qualizza
 * Purpose: 
 *
 * This file is licensed under the Eclipse Public License (EPL)
 */

#include <orthocut.hpp>
#include <CglCutGenerator.hpp>
#include <OsiSolverInterface.hpp>
#include <dsyevx_wrapper.hpp>
#include <tracer.hpp>
#include <misc_util.hpp>


void orthoCutGen(const double *sol, int n, OsiCuts &cs, double *z, double *w, int m, Tracer *tracer) {
	// Y = / 1 x \
	//     \xT X /
	// D = diagonal matrix of eigenvalues \lambda_i of Y
	// P = orthogonal matrix where columns are eigenvectors v_i of Y
	// Y D = Y P by definition of eigenvalue/eigenvectors
	// so Y = P D P^{-1} = P D P^T since P is orthogonal
	// let
	// Pn = matrix where columns are eigenvectors of v_i of Y with \lambda_i < 0, other colums are 0
	// Pp = matrix where columns are eigenvectors of v_i of Y with \lambda_i >= 0, other colums are 0
	// Y = Pn D Pn^T + Pp D Pp^T
	// we want to find the projection of Y on the SDP cone, which is Pp D Pp^T = Y - Pn D Pn^T

	Timer orthocut_timer;
	orthocut_timer.start();
	int origcuts = cs.sizeCuts();

	int np = n+1;
	int sol_card = n*(n+3)/2;
	
	double **Pn = new double*[np];
	for(int i=0;i<np;i++)
		Pn[i] = new double[np];
	double **PpDPpT = new double*[np];
	for(int i=0;i<np;i++)
		PpDPpT[i] = new double[np];
	double **PnDPnT = new double*[np];
	for(int i=0;i<np;i++)
		PnDPnT[i] = new double[np];
	double **Amat = new double*[np];
	for(int i=0;i<np;i++)
		Amat[i] = new double[np];
	
	// populating Amat
	Amat[0][0]=1;
	for(int i=0;i<n;i++) {
		Amat[0][i+1]		= sol[i];
		Amat[i+1][0]		= sol[i];
		Amat[i+1][i+1]		= sol[indexQ(i,i,n)];
		for(int j=i+1;j<n;j++){
			Amat[i+1][j+1]	= sol[indexQ(i,j,n)];
			Amat[j+1][i+1]	= sol[indexQ(i,j,n)];
		}
	}

	// computing Pn
	for(int i=0;i<np;i++)
		for(int j=0;j<np;j++)
			Pn[i][j] = 0.0;
	for(int i=0;i<np;i++) {
		if ((i<m) && (w[i] < 0)) {
			double *zbase = z + i * np;
			for(int j=0;j<np;j++)
				Pn[j][i] = zbase[j];
		} else
			break;
	}
	
	// computing PpDPpT
	for (int i=0;i<np;i++) {
		for (int j=0;j<np;j++) {
			PpDPpT[i][j] = Amat[i][j];
			PnDPnT[i][j] = 0.0;
			for (int k=0;k<np;k++) {
				if ((k<m) && (w[k] < 0)) {
					PpDPpT[i][j] -= w[k] * Pn[i][k] * Pn[j][k];
					PnDPnT[i][j] += w[k] * Pn[i][k] * Pn[j][k];
				}
			}
		}
	}

	double *point_on_sdp_cone = new double[sol_card];
	for(int i=0;i<n;i++) {
		point_on_sdp_cone[i] = PpDPpT[0][i+1];
		for (int j=i;j<n;j++)
			point_on_sdp_cone[indexQ(i,j,n)] = PpDPpT[i+1][j+1];
	}
	double point_on_sdp_cone00 = PpDPpT[0][0];

	double *PnDPnTvector = new double[sol_card];
	for (int i=0;i<n;i++) {
		PnDPnTvector[i] = PnDPnT[0][i+1];
		for (int j=i;j<n;j++)
			PnDPnTvector[indexQ(i,j,n)] = PnDPnT[i+1][j+1];
	}




	double *tangent_line_coeff = new double[sol_card];

	for(int i=0;i<n;i++) {
		tangent_line_coeff[i] = -2*PnDPnTvector[i];
		tangent_line_coeff[indexQ(i,i,n)] = -PnDPnTvector[indexQ(i,i,n)];
		for(int j=i+1;j<n;j++)
			tangent_line_coeff[indexQ(i,j,n)] = -2*PnDPnTvector[indexQ(i,j,n)];
	}

	//I want that the line passes through point_on_sdp_cone
	double rhs = PnDPnT[0][0] * point_on_sdp_cone00;
	for(int i=0;i<sol_card;i++)
 		rhs += point_on_sdp_cone[i] * tangent_line_coeff[i];
	
	OsiRowCut *cut = new OsiRowCut;
	double *coeff = new double [sol_card];
	int *ind   = new int [sol_card];
	int nz = 0;
	for(int i=0;i<sol_card;i++) {
		if (tangent_line_coeff[i]) {
			coeff[nz] = tangent_line_coeff[i];
			ind[nz] = i;
			nz++;
		}
	}
	cut -> setRow ( nz, ind, coeff );
	cut -> setLb (rhs);
//cut->print();

#ifdef TRACE_ORTHOCUT
	printf("orthocut violation-->%.5f\n",cut -> violated(sol));
#endif
	cs.insert (cut);

	

	delete [] PnDPnTvector;
	for(int i=0;i<np;i++)
		delete [] Pn[i];
	delete [] Pn;
	for(int i=0;i<np;i++)
		delete [] PpDPpT[i];
	delete [] PpDPpT;
	for(int i=0;i<np;i++)
		delete [] PnDPnT[i];
	delete [] PnDPnT;
	for(int i=0;i<np;i++)
		delete [] Amat[i];
	delete [] Amat;
	delete [] point_on_sdp_cone;
	delete [] tangent_line_coeff;
	delete [] coeff;
	delete [] ind;


	tracer->setOrthocutTime(orthocut_timer.time());
	tracer->setOrthocutTotalCuts(cs.sizeCuts() - origcuts);
}


