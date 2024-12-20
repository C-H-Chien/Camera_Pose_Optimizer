#include <algorithm>
#include <utility>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h> 
#include <stdlib.h> 

#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Core>
#include <Eigen/Dense>


class Photometry_Minimizer {

public:

    //> Constructor
    Photometry_Minimizer(double, double, double, double, int, int, Eigen::Matrix<double, 4, 4, Eigen::RowMajor>, \
                         Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>, \
                         Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>, \
                         Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>, \
                         Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>, \
                         Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>, \
                         Eigen::Matrix<bool, Eigen::Dynamic, Eigen::RowMajor>, \
                         Eigen::Matrix<bool, Eigen::Dynamic, Eigen::RowMajor>, \
                         bool );

    double Project_to_Current_Frame( const Eigen::Matrix<double, 4, 4, Eigen::RowMajor>& invPose );
    void run_solver( );

    void solve_systems_of_equations(const float lambda, Eigen::Matrix<double, 6 , Eigen::RowMajor>& pose_update);
    void update_camera_pose(Eigen::Matrix<double, 6 , Eigen::RowMajor>& pose, Eigen::Matrix<double, 6 , Eigen::RowMajor>& del);
    float interpolateVector(const Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>& toInterp, float x, float y, int w) const;
    Eigen::Matrix<double, 4, 4, Eigen::RowMajor> inversePoseEigen( Eigen::Matrix<double, 4, 4, Eigen::RowMajor> pose ) const;

    //> solution from the Levenburg-Marquardt optimization
    Eigen::Matrix<double, 4, 4, Eigen::RowMajor> relative_pose;

private:

    Eigen::Matrix<double, 4, 4, Eigen::RowMajor> se3ExpEigen( const Eigen::Matrix<double, 6 , Eigen::RowMajor>& input );
    Eigen::Matrix<double, 6 , Eigen::RowMajor> se3LogEigen( Eigen::Matrix<double, 4, 4, Eigen::RowMajor> Rt );

    double fx;
    double fy;
    double cx;
    double cy;
    int img_w;
    int img_h;
    bool mex_debug;
    
    // Image Vectors and residual vectors
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_im1;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_im2;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_im1Final;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_im2Final;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_residual;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_rsquared;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_weights;
    Eigen::Matrix<float, Eigen::Dynamic, 6, Eigen::RowMajor> m_Jacobian;

    Eigen::Matrix<float, 3, Eigen::Dynamic, Eigen::RowMajor>    m_X3D_Cur_Frame;
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>    m_X3D;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::RowMajor>        m_finalMask;          //> final edge list (just a list of all 1s)
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::RowMajor>        m_edgeMask;

    // Vectors of Image Gradients
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_gx;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_gy;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_gxFinal;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_gyFinal;

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_XFinal;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_YFinal;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_ZFinal;

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_warpedX;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_warpedY;

    //> Solver parameters
    const int MAX_ITERATIONS = 50;
    float lambda;
    const float INF_F = std::numeric_limits<float>::infinity();
    const float HUBER_THRESH = 5.f;
    const float LAMBDA_MIN = 0.f;
    const float LAMBDA_MAX = 0.2f;
    const float LAMBDA_UPDATE_FACTOR = 0.5f;
    const float MIN_GRADIENT_THRESH = 50.f;

    //> Some camera pose update parameters
    const double MIN_TRANSLATION_UPDATE = 1.e-8;
    const double MIN_ROTATION_UPDATE = 1.e-8;
    const double EPSILON = 1.e-8;
};

