#include "mex.h"
#include "photometry_minimizer.hpp"

Photometry_Minimizer::Photometry_Minimizer(double FX, double FY, double CX, double CY, \
                                           int img_width, int img_height, Eigen::Matrix<double, 4, 4, Eigen::RowMajor> relative_pose, \
                                           Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_im_current, \
                                           Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_im_reference_edges, \
                                           Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_gx, \
                                           Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> m_gy, \
                                           Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> m_X3D, \
                                           Eigen::Matrix<bool, Eigen::Dynamic, Eigen::RowMajor> m_edgeMask, \
                                           Eigen::Matrix<bool, Eigen::Dynamic, Eigen::RowMajor> m_finalMask, \
                                           bool mex_debug )
    : fx(FX), fy(FY), cx(CX), cy(CY), img_w(img_width), img_h(img_height), relative_pose(relative_pose), \
      m_im2(m_im_current), m_im1Final(m_im_reference_edges), m_gx(m_gx), m_gy(m_gy), m_X3D(m_X3D), \
      m_edgeMask(m_edgeMask), m_finalMask(m_finalMask), mex_debug(mex_debug)
{
    //> Initialize solver hyper-parameters
    lambda = 0.f;
    // mexPrintf("img_w = %d, img_h = %d\n", img_w, img_h);
}

void Photometry_Minimizer::run_solver( ) {

    float error_last = INF_F;
    float error = error_last;

    //> Levenberg-Marquardt algorithm
    for(int i = 0; i < MAX_ITERATIONS; ++i)
    {
        error_last = error;
        
        Eigen::Matrix<double, 4, 4, Eigen::RowMajor> invPose = inversePoseEigen(relative_pose);
        error = Project_to_Current_Frame( invPose );
        // mexPrintf("error = %f\n", error);

        if( error < error_last )
        {
            //> Update relative pose
            Eigen::Matrix<double, 6 , Eigen::RowMajor> del;
            solve_systems_of_equations(lambda, del);
            
            if( (del.segment<3>(0)).dot(del.segment<3>(0)) < MIN_TRANSLATION_UPDATE & 
                (del.segment<3>(3)).dot(del.segment<3>(3)) < MIN_ROTATION_UPDATE    )
                break;

            //> Convert lie algebra living in the lie group to transformation matrix
            Eigen::Matrix<double, 4, 4, Eigen::RowMajor> delMat = se3ExpEigen(del);
            relative_pose = delMat * relative_pose;

            //> Update lambda
            if(lambda <= LAMBDA_MAX)
                lambda = LAMBDA_MIN;
            else
                lambda *= LAMBDA_UPDATE_FACTOR;
        }
        else
        {
            if(lambda == LAMBDA_MIN)
                lambda = LAMBDA_MAX;
            else
                lambda *= LAMBDA_UPDATE_FACTOR;
        }
    }
    mexPrintf("error = %f\n", error);
    // std::cout << "====== Final Optimized Relative Pose =====" << std::endl;
    // std::cout << relative_pose.getPoseMatrix() << std::endl;
}

double Photometry_Minimizer::Project_to_Current_Frame( const Eigen::Matrix<double, 4, 4, Eigen::RowMajor>& invPose )
{
    //> Get the relative rotation and translation
    Eigen::Matrix<float,3,3> R = (invPose.block<3,3>(0,0)).cast<float>() ;
    Eigen::Matrix<float,3,1> t = (invPose.block<3,1>(0,3)).cast<float>() ;
    
    //> Transform 3D edge points from reference frame coordinate to the current frame coordinate
    m_X3D_Cur_Frame.resize(Eigen::NoChange, m_X3D.rows());
    m_X3D_Cur_Frame = R * m_X3D.transpose() + t.replicate(1, m_X3D.rows() );

    //> Project 3D edge point to 2D current frame
    m_warpedX.resize(m_X3D.rows());
    m_warpedY.resize(m_X3D.rows());
    m_warpedX = (fx * (m_X3D_Cur_Frame.row(0)).array() / (m_X3D_Cur_Frame.row(2)).array() ) + cx;
    m_warpedY = (fy * (m_X3D_Cur_Frame.row(1)).array() / (m_X3D_Cur_Frame.row(2)).array() ) + cy;

    //> Turning the flag for invalid 3D edges
    m_finalMask = m_edgeMask;

    m_finalMask = ( m_X3D_Cur_Frame.row(2).transpose().array() <= 0.f ).select(0, m_finalMask);
    m_finalMask = ( m_X3D.col(2).array() <= 0.f ).select(0, m_finalMask);
    m_finalMask = ( (m_X3D.col(2).array()).isFinite() ).select(m_finalMask, 0);
    m_finalMask = ( (m_X3D_Cur_Frame.row(2).transpose().array()).isFinite() ).select(m_finalMask, 0);
    
    //> Check new projected x coordinates are: 0 <= x < img_w-1
    m_finalMask = (m_warpedX.array() < 0.f).select(0, m_finalMask);
    m_finalMask = (m_warpedX.array() >= img_w-2).select(0, m_finalMask);
    //> Check new projected x coordinates are: 0 <= y < img_h-1
    m_finalMask = (m_warpedY.array() >= img_h-2).select(0, m_finalMask);
    m_finalMask = (m_warpedY.array() < 0.f).select(0, m_finalMask);

    m_finalMask = (m_warpedY.array().isFinite()).select(m_finalMask, 0);
    
    size_t numElements = (m_finalMask.array() != 0).count();
    // mexPrintf("numElements = %d\n", numElements);
    m_gxFinal.resize(numElements);
    m_gyFinal.resize(numElements);
    m_im1.resize(numElements);
    m_im2Final.resize(numElements);
    m_XFinal.resize(numElements);
    m_YFinal.resize(numElements);
    m_ZFinal.resize(numElements);

    size_t idx = 0;
    for(int i = 0; i < m_finalMask.rows(); ++i)
    {
        if(m_finalMask[i] != 0)
        {
            m_gxFinal[idx]  = interpolateVector( m_gx, m_warpedX[i], m_warpedY[i], img_w );
            m_gyFinal[idx]  = interpolateVector( m_gy, m_warpedX[i], m_warpedY[i], img_w) ;
            m_im1[idx]      = m_im1Final[i];
            m_im2Final[idx] = interpolateVector( m_im2, m_warpedX[i], m_warpedY[i], img_w );
            m_XFinal[idx]   = m_X3D_Cur_Frame(0,i);
            m_YFinal[idx]   = m_X3D_Cur_Frame(1,i);
            m_ZFinal[idx]   = m_X3D_Cur_Frame(2,i);
            ++idx;
        }
    }
    
    //> Calculate A and b matrices
    m_residual.resize(numElements);
    m_rsquared.resize(numElements);
    m_weights.resize(numElements);

    //> The residuals and the corresponding weights 
    m_residual = ( m_im1.array() - m_im2Final.array() );
    m_rsquared = m_residual.array() * m_residual.array();
    m_weights = Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>::Ones(numElements);
    m_weights = ( ( (m_residual.array()).abs() ) > HUBER_THRESH ).select( HUBER_THRESH / (m_residual.array()).abs() , m_weights);

    //> Return the entire weighted average residual
    return ( (m_weights.array() * m_rsquared.array()).sum() / (float) numElements );
}

float Photometry_Minimizer::interpolateVector(const Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor>& toInterp, float x, float y, int w) const
{
    int xi = (int) x;
	int yi = (int) y;
	float dx = x - xi;
	float dy = y - yi;
	float dxdy = dx * dy;
    int topLeft = w * yi + xi;
    int topRight = topLeft + 1;
    int bottomLeft = topLeft + w;
    int bottomRight= bottomLeft + 1;
  
    //               x                x+1
    //       ======================================
    //  y    |    topLeft      |    topRight      |
    //       ======================================
    //  y+w  |    bottomLeft   |    bottomRight   |
    //       ======================================
    return  dxdy * toInterp[bottomRight]
	        + (dy - dxdy) * toInterp[bottomLeft]
	        + (dx - dxdy) * toInterp[topRight]
			+ (1.f - dx - dy + dxdy) * toInterp[topLeft];
}

void Photometry_Minimizer::solve_systems_of_equations(const float lambda, Eigen::Matrix<double, 6 , Eigen::RowMajor>& pose_update)
{
    size_t numElements = m_im2Final.rows();
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::RowMajor> Z2 = m_ZFinal.array() * m_ZFinal.array();

    m_Jacobian.resize(numElements, Eigen::NoChange);
    m_Jacobian.col(0) =  m_weights.array() * fx * ( m_gxFinal.array() / m_ZFinal.array() );

    m_Jacobian.col(1) =  m_weights.array() * fy * ( m_gyFinal.array() / m_ZFinal.array() );

    m_Jacobian.col(2) = - m_weights.array()* ( fx * ( m_XFinal.array() * m_gxFinal.array() ) + fy * ( m_YFinal.array() * m_gyFinal.array() ) )
                        / ( Z2.array() );

    m_Jacobian.col(3) = - m_weights.array() * ( fx * m_XFinal.array() * m_YFinal.array() * m_gxFinal.array() / Z2.array()
                         + fy *( 1.f + ( m_YFinal.array() * m_YFinal.array() / Z2.array() ) ) * m_gyFinal.array() );

    m_Jacobian.col(4) = m_weights.array() * ( fx * (1.f + ( m_XFinal.array() * m_XFinal.array() / Z2.array() ) ) * m_gxFinal.array() 
                        + fy * ( m_XFinal.array() * m_YFinal.array() * m_gyFinal.array() ) / Z2.array() );

    m_Jacobian.col(5) = m_weights.array() * ( -fx * ( m_YFinal.array() * m_gxFinal.array() ) + fy * ( m_XFinal.array() * m_gyFinal.array() ) )
                        / m_ZFinal.array();
    
    m_residual.array() *= m_weights.array();
    
    pose_update = -( (m_Jacobian.transpose() * m_Jacobian).cast<double>() ).ldlt().solve( (m_Jacobian.transpose() * m_residual).cast<double>() );   
}

Eigen::Matrix<double, 4, 4, Eigen::RowMajor> Photometry_Minimizer::inversePoseEigen( Eigen::Matrix<double, 4, 4, Eigen::RowMajor> pose ) const
{
    Eigen::Matrix<double, 4, 4, Eigen::RowMajor> invPose = pose;
    invPose.block<3,3>(0,0).transposeInPlace();
    invPose.block<3,1>(0,3) = -invPose.block<3,3>(0,0) * invPose.block<3,1>(0,3);
    return invPose;
}

Eigen::Matrix<double, 4, 4, Eigen::RowMajor> Photometry_Minimizer::se3ExpEigen( const Eigen::Matrix<double, 6 , Eigen::RowMajor>& input )
{
    //> Size of input must be [6 x 1] or [1 x 6]
    Eigen::MatrixXd lie(4, 4);
    // Eigen::Matrix<double, 4, 4,Eigen::RowMajor> lie;
    lie << 0.,           -input[5],  input[4],   input[0],
           input[5],      0.,        -input[3],  input[1],
           -input[4],     input[3],   0.,        input[2],
           0.,           0.,        0.,        0.;

    //> Exponentiate
    Eigen::Matrix<double, 4, 4, Eigen::RowMajor> Rt( lie.exp() );
    return Rt;
}

Eigen::Matrix<double, 6 , Eigen::RowMajor> Photometry_Minimizer::se3LogEigen( Eigen::Matrix<double, 4, 4, Eigen::RowMajor> Rt )
{
    //> Size of input must be [6 x 1] or [1 x 6]
    Eigen::Matrix<double, 6 , Eigen::RowMajor> logMap;
    Eigen::Matrix<double,3,3, Eigen::RowMajor> R = Rt.block<3,3>(0,0);
    Eigen::Matrix<double,3, Eigen::RowMajor> t = Rt.block<3,1>(0,3);
    Eigen::Matrix<double,3,3, Eigen::RowMajor> A = 0.5*(R - R.transpose());
    double s = sqrt( A(2,1)*A(2,1) + A(0,2) * A(0,2) + A(1,0) * A(1,0) );
    double c = 0.5*(R.trace() - 1);
    if(s < EPSILON & (c-1.) < EPSILON)
    {
        logMap << t , 0. , 0., 0.;
        return logMap;
    }
    else
    {
        double theta = atan2(s,c);
        Eigen::Matrix<double,3, Eigen::RowMajor> r;
        r << A(2,1) , A(0,2) , A(1,0);
        r = (theta/s)*r;
        Eigen::Matrix<double,3,3, Eigen::RowMajor> V;
        Eigen::Matrix<double,3,3, Eigen::RowMajor> wx;
        Eigen::Matrix<double,3, Eigen::RowMajor> tp;
        wx << 0. , -r[2], r[1],
              r[2], 0,    -r[0],
              -r[1], r[0], 0.;
        V.setIdentity();
        V = V + ((1.-c)/ pow(theta,2.)) * wx  + ((theta- s)/pow(theta,3.))*wx*wx;
        tp = V.ldlt().solve(t);
        logMap.head<3>() = tp;
        logMap.tail<3>() = r;
        return logMap;
    }
}
