#include <algorithm>
#include <utility>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h> 
#include <stdlib.h> 

#include <Eigen/Core>
#include <Eigen/Dense>
#include "mex.h"

#include "photometry_minimizer.hpp"

//> main mex function
void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] ) {

    //> Get data from MATLAB
    double* current_frame               = (double*) mxGetData(pr[0]);
    double* reference_frame             = (double*) mxGetData(pr[1]);
    double* reference_edge_mask         = (double*) mxGetData(pr[2]);
    double* reference_edge_3D_location  = (double*) mxGetData(pr[3]);
    double* reference_frame_gx          = (double*) mxGetData(pr[4]);
    double* reference_frame_gy          = (double*) mxGetData(pr[5]);
    double* intrinsic_matrix            = (double*) mxGetData(pr[6]);
    double* init_relative_pose          = (double*) mxGetData(pr[7]);
    int     num_of_edges                = (int) mxGetScalar(pr[8]);
    mexPrintf("Number of edges = %d\n", num_of_edges);
    int     img_width                   = (int) mxGetScalar(pr[9]);
    int     img_height                  = (int) mxGetScalar(pr[10]);
    bool    mexFunction_debug           = (bool) mxGetScalar(pr[11]);

    if (mexFunction_debug) {
        mexPrintf("Number of edges = %d\n", num_of_edges);
    }

    double fx = intrinsic_matrix[0];
    double fy = intrinsic_matrix[4];
    double cx = intrinsic_matrix[2];
    double cy = intrinsic_matrix[5];
    int full_img_size = img_width * img_height;

    //> A 4x4 initial guess of a transformation matrix 
    Eigen::Matrix<double, 4, 4, Eigen::RowMajor> relative_pose;
    relative_pose << init_relative_pose[0], init_relative_pose[1], init_relative_pose[2], init_relative_pose[3], \
                     init_relative_pose[4], init_relative_pose[5], init_relative_pose[6], init_relative_pose[7], \
                     init_relative_pose[8], init_relative_pose[9], init_relative_pose[10], init_relative_pose[11], \
                     0, 0, 0, 1;
    if (mexFunction_debug) {
        mexPrintf("fx = %f, fy = %f, cx = %f, cy = %f\n", fx, fy, cx, cy);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                mexPrintf("%f\t", relative_pose(i*4 + j));
            }
            mexPrintf("\n");
        }
        for (int i = num_of_edges-2; i > num_of_edges-10; i--)
            mexPrintf("(%f, %f, %f)\n", reference_edge_3D_location[i*3], reference_edge_3D_location[i*3+1], reference_edge_3D_location[i*3+2]);
    }

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>       m_im_reference;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>       m_im_current;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>       m_im_reference_edges; //> m_im1Final;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>       m_gx;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>       m_gy;
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>    m_X3D;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::RowMajor>        m_finalMask;          //> final edge list (just a list of all 1s)
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::RowMajor>        m_edgeMask;           //> fedge list (just a list of all 1s)

    m_im_reference.resize(full_img_size);
    m_im_current.resize(full_img_size);
    m_gx.resize(full_img_size);
    m_gy.resize(full_img_size);
    m_im_reference_edges.resize(num_of_edges);
    m_edgeMask.resize(num_of_edges);
    m_finalMask.resize(num_of_edges);
    m_X3D.resize(num_of_edges ,Eigen::NoChange);

    //> Assigning input values to the Eigen matrices of *images*
    int idx = 0;
    int test_count = 0;
    float edge_X, edge_Y, edge_Z;
    for (int i = 0; i < full_img_size; i++) {
        m_im_current[i]     = current_frame[i];
        m_im_reference[i]   = reference_frame[i];
        m_gx[i]             = reference_frame_gx[i];
        m_gy[i]             = reference_frame_gy[i];

        //> if it's an edge, make arrays specifically for edges
        test_count += reference_edge_mask[i];
        if (reference_edge_mask[i] == 1) {
            edge_X          = reference_edge_3D_location[idx*3 + 0];
            edge_Y          = reference_edge_3D_location[idx*3 + 1];
            edge_Z          = reference_edge_3D_location[idx*3 + 2];
            m_X3D.row(idx)  << edge_X, edge_Y, edge_Z;
            m_im_reference_edges[idx] = reference_frame[i];
            m_edgeMask[idx] = (bool)(reference_edge_mask[i]);
            m_finalMask[idx] = (bool)(reference_edge_mask[i]);
            idx++;
        }
    }
    if (mexFunction_debug) {
        mexPrintf("test_count = %d\n", test_count);
        mexPrintf("idx = %d\n", idx);
        for (int i = idx-2; i > idx-10; i--)
            mexPrintf("(%f, %f, %f)\n", (m_X3D.row(i))(0), (m_X3D.row(i))(1), (m_X3D.row(i))(2));
    }

    Photometry_Minimizer Camera_Pose_Solver(fx, fy, cx, cy, img_width, img_height, relative_pose, \
                                            m_im_current, m_im_reference_edges, m_gx, m_gy, m_X3D, m_edgeMask, m_finalMask, mexFunction_debug);
    Camera_Pose_Solver.run_solver( );
    
    //> Optimized relative pose
    if (mexFunction_debug) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                mexPrintf("%f\t", Camera_Pose_Solver.relative_pose(i*4 + j));
            }
            mexPrintf("\n");
        }
    }
    
    //> Output reprojection errors to MATLAB
    if ( nl > 0 ) {
        const mwSize optimized_pose_dim[1] = { mwSize(16) };
        pl[0] = mxCreateNumericArray(1, optimized_pose_dim, mxDOUBLE_CLASS, mxREAL);
        double* output_optimized_pose = (double*) mxGetData(pl[0]);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                output_optimized_pose[i*4 + j] = Camera_Pose_Solver.relative_pose(i*4 + j);
            }
        }
    }
    // camera_pose.updateKeyFramePose(keyframe_pose.getPoseMatrix(), relative_pose.getPoseMatrix());
}
