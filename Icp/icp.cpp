/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

Authors: Andreas Geiger

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

#include "../include/ICP/icp.h"

using namespace std;

Icp::Icp (float *M,const int32_t M_num,const int32_t dim) :
  dim(dim), max_iter(10), min_delta(1e-4),indis(-1.0),template_num(0) {


  // check for correct dimensionality
  if (dim!=2 && dim!=3) {
    cout << "ERROR: LIBICP works only for data of dimensionality 2 or 3" << endl;
    M_tree = 0;
    return;
  }

  // check for minimum number of points
  if (M_num<5) {
    cout << "ERROR: LIBICP works only with at least 5 model points" << endl;
    M_tree = 0;
    return;
  }

  // copy model points to M_data
  M_data.resize(boost::extents[M_num][dim]);
  for (int32_t m=0; m<M_num; m++)
    for (int32_t n=0; n<dim; n++)
      M_data[m][n] = M[m*dim+n];

  // build a kd tree from the model point cloud
  M_tree = new kdtree::KDTree(M_data);
}

Icp::~Icp () {
  if (M_tree)
    delete M_tree;
}

void Icp::fit (float *T,const int32_t T_num,Matrix &R,Matrix &t,const double indist) {

  // make sure we have a model tree
  if (!M_tree) {
    cout << "ERROR: No model available." << endl;
    return;
  }

  // check for minimum number of points
  if (T_num<5) {
    cout << "ERROR: Icp works only with at least 5 template points" << endl;
    return;
  }


  // set active points
 template_num = T_num;//template data number
 indis = indist;
//save template data
  T_data = Matrix(3,T_num);
  for(int i = 0 ; i < T_num; i++)
  {
      T_data.val[0][i] = T[i*3+0];
      T_data.val[1][i] = T[i*3+1];
      T_data.val[2][i] = T[i*3+2];

  }
  ///applytion R|T

 this->T_data = R * this->T_data;
 for(int i= 0; i < template_num; i++)
 {
    T_data.val[0][i] = T_data.val[0][i] + t.val[0][0];
    T_data.val[1][i] = T_data.val[1][i] + t.val[1][0];
    T_data.val[2][i] = T_data.val[2][i] + t.val[2][0];

 }
  vector<int32_t> active;
  if (indist<=0) {
    active.clear();
    for (int32_t i=0; i<T_num; i++)
      active.push_back(i);
  } else {
    active = getInliers(T,T_num,R,t,indist);
  }

  // run icp

  fitIterate(T,T_num,R,t,active);

}

void Icp::fitIterate(float *T,const int32_t T_num,Matrix &R,Matrix &t,const std::vector<int32_t> &active) {

  // check if we have at least 5 active points
  if (active.size()<5)
    return;

  // iterate until convergence
  for (int32_t iter=0; iter<max_iter; iter++)
    if (fitStep(T,T_num,R,t,active)<min_delta)
      break;
}
