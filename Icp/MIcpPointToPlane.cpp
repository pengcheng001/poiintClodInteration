/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

Authors: Andreas Geiger

openMP support by Manolis Lourakis, Foundation for Research & Technology - Hellas, Heraklion, Greece

libicp is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or any later version.

libicp is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libicp; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA
*/

#ifdef _OPENMP
#include <omp.h>
#endif
#include<fstream>

#include"../include/ICP/MIcpPointToPlane.h"
using namespace std;

// Also see (3d part): "Least-Squares Fitting of Two 3-D Point Sets" (Arun, Huang and Blostein)
double MIcpPointToPlane::fitStep (float *T,const int32_t T_num,Matrix &R,Matrix &t,const std::vector<int32_t> &actives) {

  int i;
  vector<int32_t> active;
      if (indis<=0) {
        active.clear();
        for (int32_t i=0; i<T_num; i++)
          active.push_back(i);
      } else {
        active = getInliers(T,T_num,R,t,indis);
      }
  int nact = (int)active.size();

  // init matrix for point correspondences
  Matrix p_m(dim,nact); // model
  Matrix p_t(dim,nact);
  Matrix resids(nact,1);//
  std::vector<double> resids_t;


    // extract matrix and translation vector
    /*
    double r00 = R.val[0][0]; double r01 = R.val[0][1]; double r02 = R.val[0][2];
    double r10 = R.val[1][0]; double r11 = R.val[1][1]; double r12 = R.val[1][2];
    double r20 = R.val[2][0]; double r21 = R.val[2][1]; double r22 = R.val[2][2];
    double t0  = t.val[0][0]; double t1  = t.val[1][0]; double t2  = t.val[2][0];*/
    // establish correspondences
  for (i=0; i<nact; i++) {
      // kd tree query + result
      std::vector<float> query(dim);
      kdtree::KDTreeResultVector result;

      // get index of active point
      int32_t idx = active[i];

      // transform point according to R|t
      double sx = T_data.val[0][idx] ; query[0] = (float)sx;
      double sy = T_data.val[1][idx] ; query[1] = (float)sy;
      double sz = T_data.val[2][idx] ; query[2] = (float)sz;

	  
      // search nearest neighbor
      M_tree->n_nearest(query,1,result);
      // set model point
      p_m.val[0][i] = M_tree->the_data[result[0].idx][0];
      p_m.val[1][i] = M_tree->the_data[result[0].idx][1];
      p_m.val[2][i] = M_tree->the_data[result[0].idx][2];

      // set template point
      p_t.val[0][i] = query[0];
      p_t.val[1][i] = query[1];
      p_t.val[2][i] = query[2];

	 
      resids.val[i][0] = sqrt(
                         pow((p_m.val[0][i] - p_t.val[0][i]),2.0)+
                         pow((p_m.val[1][i] - p_t.val[1][i]),2.0)+
                         pow((p_m.val[2][i] - p_t.val[2][i]),2.0)
                         );
     resids_t.push_back(resids.val[i][0]);
	/* std::cout<<"model:"<<p_m.val[0][i]<<", "<<p_m.val[1][i]<<", "<<p_m.val[2][i]<<std::endl;
	 std::cout<<"template:"<<p_t.val[0][i]<<", "<<p_t.val[1][i]<<", "<<p_t.val[2][i]<<", "<<std::endl;
	 std::cout<<"res"<<resids.val[i][0]<<std::endl;*/
    }
     std::sort(resids_t.begin(), resids_t.end(),std::greater<double>());
    ///solve
      double median = 0;
      double max_res = resids_t[0];
   // std::vector<double> wghs;
    if(resids_t.size()== 0) /// resids.size() / 0.0 == 0
    {
        median  = resids_t[resids_t.size() / 2];
    }else
    {
        median = resids_t[(resids_t.size())/2];
    }
    double kRob = 1.9*median;
    if(kRob < 1e-6*max_res)
        kRob = 0.3*max_res;
    else if(abs(max_res - kRob) < 1e-13)
        kRob = 1;
   unsigned int cont = 0;
   long  double mean = 0.0;
    for(int i = 0 ; i < nact; i++)
    {
        if(resids_t[i] < 1.5*kRob){
		    mean += pow(resids_t[i],2.0);
			cont++;
		}
    }
    mean /= cont;
    /////////////////////////////
    Matrix wghs(nact,1);
    Matrix med;
    Matrix mem;
	method = 1;
    switch(method)
    {
    case 0: ///Huber
        kRob = 2.0138 * kRob;
        for(int i = 0; i <nact; i++)
        {
            if(resids.val[i][0] < kRob)
                wghs.val[i][0] = 1;
            else
                wghs.val[i][0] = kRob/resids.val[i][0];

        }
		break;
    case 1:///Tukey's bi-weight
        kRob = 7.0589 * kRob;
        for(int i = 0; i < nact; i++)
        {
            if(resids.val[i][0] < kRob)
                wghs.val[i][0] = pow(1.0 -  pow(resids.val[i][0]/kRob,2.0),2.0);
            else
                wghs.val[i][0] = 0.0;
			
        }
	
		break;
    case 2:///% Cauchy
        kRob = 4.3040*kRob;
        for(int i = 0; i < nact; i++)
        {
            wghs.val[i][0] = 1.0 /( 1 +  pow(resids.val[i][0]/kRob, 2.0));
        }
		break;
    case 3: ////Welsch
         kRob = 4.7536*kRob;
        for(int i = 0; i < nact; i++)
        {
            wghs.val[i][0] = exp( 0 -  pow(resids.val[i][0]/kRob, 2.0));
        }
		break;
    default:////
       kRob = 2.0138 * kRob;
        for(int i = 0; i < nact; i++)
        {
            if(resids.val[i][0] < kRob)
                wghs.val[i][0] = 1;
            else
                wghs.val[i][0] = kRob/resids.val[i][0];
        }
		break;

    }
    double sumWghs = 0.0;

    for(int i = 0; i < nact; i++)
    {
        sumWghs += wghs.val[i][0];
    }
	
    med = (p_t*wghs)/sumWghs;
    mem = (p_m*wghs)/sumWghs;


    for(int i = 0; i < nact; i++)
    {
        p_t.val[0][i] = p_t.val[0][i]*wghs.val[i][0];
        p_t.val[1][i] = p_t.val[1][i]*wghs.val[i][0];
        p_t.val[2][i] = p_t.val[2][i]*wghs.val[i][0];

    }
    Matrix C = p_t*(~p_m) - (med*sumWghs)*(~mem);
    Matrix U,W,V,Ri,Ti;
    C.svd(U,W,V);
    Ri = V*(~U);
    if(Ri.det() < 0)
    {
        V.val[0][2] = - V.val[0][2];
        V.val[1][2] = - V.val[1][2];
        V.val[2][2] = - V.val[2][2];
        Ri = V*(~U);
    }

    Ti = mem - Ri*med;
    this->T_data = Ri* this->T_data;
	
    for(int i= 0; i < template_num; i++)
    {
        T_data.val[0][i] = T_data.val[0][i] + Ti.val[0][0];
        T_data.val[1][i] = T_data.val[1][i] + Ti.val[1][0];
        T_data.val[2][i] = T_data.val[2][i] + Ti.val[2][0];

    }
    R = Ri*R;
    t = Ri*t + Ti;
    return mean;

}

std::vector<int32_t> MIcpPointToPlane::getInliers (float *T,const int32_t T_num,const Matrix &R,const Matrix &t,const double indist) {

   // init inlier vector + query point + query result
  vector<int32_t>            inliers;
  std::vector<float>         query(dim);
  kdtree::KDTreeResultVector neighbor;

  // dimensionality 2
   {


    // check for all points if they are inliers
    for (int32_t i=0; i<T_num; i++) {

      // transform point according to R|t
      double sx = T_data.val[0][i] + T_data.val[1][i] + T_data.val[2][i] ; query[0] = (float)sx;
      double sy = T_data.val[1][i] + T_data.val[1][i] + T_data.val[2][i] ; query[1] = (float)sy;
      double sz = T_data.val[2][i] + T_data.val[1][i] + T_data.val[1][i] ; query[2] = (float)sz;
      // search nearest neighbor
      M_tree->n_nearest(query,1,neighbor);

      // model point
      double dx = M_tree->the_data[neighbor[0].idx][0];
      double dy = M_tree->the_data[neighbor[0].idx][1];
      double dz = M_tree->the_data[neighbor[0].idx][2];

      // model point normal
      double nx = M_normal[neighbor[0].idx*3+0];
      double ny = M_normal[neighbor[0].idx*3+1];
      double nz = M_normal[neighbor[0].idx*3+2];

      // check if it is an inlier
      if ((sx-dx)*nx+(sy-dy)*ny+(sz-dz)*nz<indist)
        inliers.push_back(i);
    }
  }

  // return vector with inliers
  return inliers;
}

void MIcpPointToPlane::computeNormal (const kdtree::KDTreeResultVector &neighbors,double *M_normal,const double flatness) {

  // dimensionality 2
  if (dim==2) {

    // extract neighbors
    Matrix P(neighbors.size(),2);
    Matrix mu(1,2);
    for (uint32_t i=0; i<neighbors.size(); i++) {
      double x = M_tree->the_data[neighbors[i].idx][0];
      double y = M_tree->the_data[neighbors[i].idx][1];
      P.val[i][0] = x;
      P.val[i][1] = y;
      mu.val[0][0] += x;
      mu.val[0][1] += y;
    }

    // zero mean
    mu       = mu/(double)neighbors.size();
    Matrix Q = P - Matrix::ones(neighbors.size(),1)*mu;

    // principal component analysis
    Matrix H = ~Q*Q;
    Matrix U,W,V;
    H.svd(U,W,V);

    // normal
    M_normal[0] = U.val[0][1];
    M_normal[1] = U.val[1][1];

  // dimensionality 3
  } else {

    // extract neighbors
    Matrix P(neighbors.size(),3);
    Matrix mu(1,3);
    for (uint32_t i=0; i<neighbors.size(); i++) {
      double x = M_tree->the_data[neighbors[i].idx][0];
      double y = M_tree->the_data[neighbors[i].idx][1];
      double z = M_tree->the_data[neighbors[i].idx][2];
      P.val[i][0] = x;
      P.val[i][1] = y;
      P.val[i][2] = z;
      mu.val[0][0] += x;
      mu.val[0][1] += y;
      mu.val[0][2] += z;
    }

    // zero mean
    mu       = mu/(double)neighbors.size();
    Matrix Q = P - Matrix::ones(neighbors.size(),1)*mu;

    // principal component analysis
    Matrix H = ~Q*Q;
    Matrix U,W,V;
    H.svd(U,W,V);

    // normal
    M_normal[0] = U.val[0][2];
    M_normal[1] = U.val[1][2];
    M_normal[2] = U.val[2][2];
  }
}

double* MIcpPointToPlane::computeNormals (const int32_t num_neighbors,const double flatness) {
  double *M_normal = (double*)malloc(M_tree->N*dim*sizeof(double));
  kdtree::KDTreeResultVector neighbors;
  for (int32_t i=0; i<M_tree->N; i++) {
    M_tree->n_nearest_around_point(i,0,num_neighbors,neighbors);
    if (dim==2) computeNormal(neighbors,M_normal+i*2,flatness);
    else        computeNormal(neighbors,M_normal+i*3,flatness);
  }
  return M_normal;
}

