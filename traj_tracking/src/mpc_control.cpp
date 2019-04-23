#include<iostream>
#include"mpc_control.h"
#include "Eigen/LU"


using namespace std;
Eigen::MatrixXd powMatrix(Eigen::MatrixXd ma,int n)
{
    Eigen::MatrixXd result_matrix=Eigen::MatrixXd::Identity(ma.rows(),ma.cols());
    for(int i=0;i<n;i++)
    {
        result_matrix*=ma;
    }
    return result_matrix;
}

MPCControl::MPCControl(double sample_time,double car_length,int Np,int Nc)
{
    //控制时域和预测时域
    Np_=Np;
    Nc_=Nc;
    sample_time_=sample_time;
    //状态变量个数和控制变量个数
    basic_state_size_=3;
    basic_control_size_=2;        
    //松弛因子
    relax_factor_=2.0;
    relax_factor_max_=4.0;
    //汽车参数
    car_length_=car_length;

    //约束参数
    //转角单个方向的范围
    wheel_single_direction_max_degree_ = 130.0/180.0*M_PI;
    //最大转角转速
    wheel_single_direction_max_vel_ = 130.0/180.0*M_PI;
    //最大速度
    vx_max_ = 300;
    //最大纵向加速度
    ax_max_ = 200;
    //当前位置误差
    vector_error_=Eigen::MatrixXd::Zero(basic_state_size_+basic_control_size_,1);

    // parameters for mpc solver; number of iterations
    mpc_max_iteration_ = 500;
    // parameters for mpc solver; threshold for computation
    mpc_eps_ = 0.05;

    real_control_vel_.push_back(0.0);
    real_control_delta_.push_back(0.0);

}


MPCControl:: ~ MPCControl()
{
    ;
}

void MPCControl::updateMatrix(double vx_ref,double phi_ref,double delta_ref)
{
    //更新线性化的vehicle状态空间方程，与离散方程
    matrix_a_=Eigen::MatrixXd::Zero(basic_state_size_,basic_state_size_);
    matrix_a_(0,2)=-vx_ref*sin(phi_ref);
    matrix_a_(1,2)=vx_ref*sin(phi_ref);

    matrix_ad_=Eigen::MatrixXd::Zero(basic_state_size_,basic_state_size_);
    matrix_ad_(0,2)=-vx_ref*sin(phi_ref)*sample_time_;
    matrix_ad_(1,2)=vx_ref*cos(phi_ref)*sample_time_;
    matrix_ad_(0,0)=1.0;
    matrix_ad_(1,1)=1.0;
    matrix_ad_(2,2)=1.0;


    matrix_b_=Eigen::MatrixXd::Zero(basic_state_size_,basic_control_size_);
    matrix_b_(0,0)=cos(phi_ref);
    matrix_b_(1,0)=sin(phi_ref);
    matrix_b_(2,0)=tan(delta_ref)/car_length_;
    matrix_b_(2,1)=vx_ref/(car_length_*cos(delta_ref)*cos(delta_ref));

    matrix_bd_=matrix_b_*sample_time_;

    matrix_c_=Eigen::MatrixXd::Identity(basic_state_size_,basic_state_size_);

    //更新状态空间方程（考虑将控制量也作为状态变量）
    matrix_A_.resize(basic_control_size_+basic_state_size_,basic_control_size_+basic_state_size_);
    matrix_A_<<matrix_ad_,matrix_bd_,Eigen::MatrixXd::Zero(basic_control_size_,basic_state_size_),Eigen::MatrixXd::Identity(basic_control_size_,basic_control_size_);

    matrix_B_.resize(basic_control_size_+basic_state_size_,basic_control_size_);
    matrix_B_<<matrix_bd_,Eigen::MatrixXd::Identity(basic_control_size_,basic_control_size_);;

    matrix_C_.resize(basic_state_size_,basic_control_size_+basic_state_size_);
    matrix_C_<<Eigen::MatrixXd::Identity(basic_state_size_,basic_state_size_),Eigen::MatrixXd::Zero(basic_state_size_,basic_control_size_);


    //cout<<matrix_A_<<endl<<endl;
    //cout<<matrix_B_<<endl<<endl;
    //cout<<matrix_C_<<endl<<endl;

    //更新预测矩阵
    matrix_Ap_=Eigen::MatrixXd::Zero(Np_*matrix_C_.rows(),matrix_A_.cols());
    matrix_Bp_=Eigen::MatrixXd::Zero(Np_*matrix_C_.rows(),Nc_*matrix_B_.cols());
    for(int i=0;i<Np_;i++)
        for(int j=0;j<Nc_;j++)
        {
            if(j==0)
            matrix_Ap_.block(i*matrix_C_.rows(),j*matrix_A_.cols(),matrix_C_.rows(),matrix_A_.cols())=matrix_C_*powMatrix(matrix_A_,i+1);

            if(i>=j)
            matrix_Bp_.block(i*matrix_C_.rows(),j*matrix_B_.cols(),matrix_C_.rows(),matrix_B_.cols())=matrix_C_*powMatrix(matrix_A_,i-j)*matrix_B_;


        }
        //cout<<matrix_C_;
    //更新权重矩阵
    matrix_q_=100*Eigen::MatrixXd::Identity(basic_state_size_*Np_,basic_state_size_*Np_);
    matrix_r_=5*Eigen::MatrixXd::Identity(basic_control_size_*Nc_,basic_control_size_*Nc_);
    //cout<<matrix_q_<<endl<<matrix_r_;

    //更新时域的参考的轨迹（参考轨迹-参考点）
    matrix_errors_ref_traj_=Eigen::MatrixXd::Zero(basic_state_size_*Np_,1);
    for(int i=0;i<Np_;i++)
    {
        matrix_errors_ref_traj_(i*basic_state_size_,0)=errors_ref_traj_[i].x;
        matrix_errors_ref_traj_(i*basic_state_size_+1,0)=errors_ref_traj_[i].y;
        matrix_errors_ref_traj_(i*basic_state_size_+2,0)=errors_ref_traj_[i].phi;
    }

    //更新二次规划矩阵
    matrix_H_=Eigen::MatrixXd::Zero(basic_control_size_*Nc_+1,basic_control_size_*Nc_+1);
    matrix_H_.block(0,0,matrix_r_.rows(),matrix_r_.cols())=(matrix_Bp_.transpose())*matrix_q_*matrix_Bp_+matrix_r_;
    matrix_H_(matrix_H_.rows()-1,matrix_H_.cols()-1)=relax_factor_;
    matrix_G_=Eigen::MatrixXd::Zero(1,basic_control_size_*Nc_+1);
    matrix_G_.leftCols(basic_control_size_*Nc_)=2*(matrix_Ap_*vector_error_).transpose()*matrix_q_*matrix_Bp_;
    //matrix_G_.leftCols(basic_control_size_*Nc_)=2*(matrix_Ap_*vector_error_-matrix_errors_ref_traj_).transpose()*matrix_q_*matrix_Bp_;
    //cout<<matrix_H_<<endl;
    //cout<<matrix_G_<<endl;


    //----------更新约束矩阵----------------------
    //更新控制增量及松弛因子的上下界
    matrix_du_ub_=Eigen::MatrixXd::Zero(basic_control_size_*Nc_+1,1);
    matrix_du_lb_=Eigen::MatrixXd::Zero(basic_control_size_*Nc_+1,1);
    for(int i=0;i<Nc_;i++)
    {
        matrix_du_ub_(basic_control_size_*i,0)=ax_max_*sample_time_;
        matrix_du_ub_(basic_control_size_*i+1,0)=wheel_single_direction_max_vel_*sample_time_;
    }
    
    matrix_du_ub_(basic_control_size_*Nc_,0)=relax_factor_max_;
    matrix_du_lb_=-1.0*matrix_du_ub_;
    //cout<<matrix_du_lb_;

    //更新控制量的上下界
    matrix_u_ub_.resize(basic_control_size_*Nc_,1);
    matrix_u_lb_.resize(basic_control_size_*Nc_,1);
    for(int i=0;i<Nc_;i++)
    {
        matrix_u_ub_(basic_control_size_*i,0)=vx_max_*sample_time_;
        matrix_u_ub_(basic_control_size_*i+1,0)=wheel_single_direction_max_degree_*sample_time_;
    }
    matrix_u_lb_=-1.0*matrix_u_ub_;
    //cout<<matrix_u_lb_;

    //更新增量的线性不等式约束
    matrix_constraint_A_=Eigen::MatrixXd::Zero(basic_control_size_*Nc_,basic_control_size_*Nc_+1);
    for(int i=0;i<Nc_;i++)
        for(int j=0;j<=i;j++)
        {
            matrix_constraint_A_.block(i*basic_control_size_,j*basic_control_size_,basic_control_size_,basic_control_size_)=Eigen::MatrixXd::Identity(basic_control_size_,basic_control_size_);
        }

    //cout<<matrix_constraint_A_;

    //更新上个时刻的控制量，以便约束矩阵A的上下界的构造
    Eigen::MatrixXd matrix_previous_control;
    matrix_previous_control.resize(basic_control_size_*Nc_,1);
    for(int i=0;i<Nc_;i++)
    {
        matrix_previous_control(basic_control_size_*i,0)=real_control_vel_.back();
        matrix_previous_control(basic_control_size_*i+1,0)=real_control_delta_.back();
    }    
    //cout<<matrix_previous_control;

    //更新约束矩阵A对应的上下界
    matrix_constraint_ubA_.resize(basic_control_size_*Nc_,1);
    matrix_constraint_lbA_.resize(basic_control_size_*Nc_,1);

    matrix_constraint_ubA_=matrix_u_ub_-matrix_previous_control;
    matrix_constraint_lbA_=matrix_u_lb_-matrix_previous_control;

    //cout<<matrix_constraint_ubA_;


}

//计算当前位置和对应轨迹参考点的差note:控制量误差也必须计算
void MPCControl::computeErrors(traj & traj_ref)
{
    vector_error_=Eigen::MatrixXd::Zero(basic_state_size_+basic_control_size_,1);
    vector_error_(0,0)=current_state_.x-traj_ref.x;
    vector_error_(1,0)=current_state_.y-traj_ref.y;
    vector_error_(2,0)=current_state_.phi-traj_ref.phi;
    vector_error_(3,0)=current_state_.vel-traj_ref.vel;
    vector_error_(4,0)=current_state_.delta-traj_ref.delta;
}

//设置机器人当前位置信息
void MPCControl::updateState(traj & current_pos)
{
    current_state_=current_pos;
}

bool MPCControl::qpSlover()
{
    USING_NAMESPACE_QPOASES
    //需要先将矩阵转换为一维数组，然后才能利用qpOASES进行求解
    double h_matrix[matrix_H_.rows()*matrix_H_.cols()];
    for(int i=0;i<matrix_H_.rows();i++)
        for(int j=0;j<matrix_H_.cols();j++)
        {
            h_matrix[i*matrix_H_.cols()+j]=matrix_H_(i,j);
            //cout<<h_matrix[i*matrix_H_.cols()+j]<<" ";
        }

    double g_matrix[matrix_G_.rows()*matrix_G_.cols()];
    for(int i=0;i<matrix_G_.rows();i++)
        for(int j=0;j<matrix_G_.cols();j++)
        {
            g_matrix[i*matrix_G_.cols()+j]=matrix_G_(i,j);
            //cout<<g_matrix[i*matrix_G_.cols()+j]<<" ";
        }

    double  lower_bound[matrix_du_lb_.rows()*matrix_du_lb_.cols()];
    for(int i=0;i<matrix_du_lb_.rows();i++)
        for(int j=0;j<matrix_du_lb_.cols();j++)
        {
            lower_bound[i*matrix_du_lb_.cols()+j]=matrix_du_lb_(i,j);
            //cout<<lower_bound[i*matrix_du_lb_.cols()+j]<<" ";
        }    

    double  upper_bound[matrix_du_ub_.rows()*matrix_du_ub_.cols()];
    for(int i=0;i<matrix_du_ub_.rows();i++)
        for(int j=0;j<matrix_du_ub_.cols();j++)
        {
            upper_bound[i*matrix_du_ub_.cols()+j]=matrix_du_ub_(i,j);
            //cout<<upper_bound[i*matrix_du_ub_.cols()+j]<<" ";
        }    

    double affine_constraint_matrix[matrix_constraint_A_.rows()*matrix_constraint_A_.cols()];
    for(int i=0;i<matrix_constraint_A_.rows();i++)
        for(int j=0;j<matrix_constraint_A_.cols();j++)
        {
            affine_constraint_matrix[i*matrix_constraint_A_.cols()+j]=matrix_constraint_A_(i,j);
            //cout<<affine_constraint_matrix[i*matrix_constraint_A_.cols()+j]<<" ";
        }      


    double  constraint_lower_bound[matrix_constraint_lbA_.rows()*matrix_constraint_lbA_.cols()];
    for(int i=0;i<matrix_constraint_lbA_.rows();i++)
        for(int j=0;j<matrix_constraint_lbA_.cols();j++)
        {
            constraint_lower_bound[i*matrix_constraint_lbA_.cols()+j]=matrix_constraint_lbA_(i,j);
            //cout<<constraint_lower_bound[i*matrix_constraint_lbA_.cols()+j]<<" ";
        }    

    double  constraint_upper_bound[matrix_constraint_ubA_.rows()*matrix_constraint_ubA_.cols()];
    for(int i=0;i<matrix_constraint_ubA_.rows();i++)
        for(int j=0;j<matrix_constraint_ubA_.cols();j++)
        {
            constraint_upper_bound[i*matrix_constraint_ubA_.cols()+j]=matrix_constraint_ubA_(i,j);
            //cout<<constraint_upper_bound[i*matrix_constraint_ubA_.cols()+j]<<" ";
        } 


    //qp求解过程
    QProblem qp_problem(basic_control_size_*Nc_+1,basic_control_size_*Nc_);
    Options my_options;
    qp_problem.setOptions(my_options);
    int nWSR=mpc_max_iteration_;//备注，一定要单独给定，不能直接把mpc_max_iteration_给init函数，否则qp不可解
    auto return_result=qp_problem.init(h_matrix,g_matrix,affine_constraint_matrix,lower_bound,upper_bound,constraint_lower_bound,constraint_upper_bound,nWSR);
    if (return_result != qpOASES::SUCCESSFUL_RETURN) {
        if (return_result == qpOASES::RET_MAX_NWSR_REACHED) {
        cout << "error :qpOASES solver failed due to reached max iteration"<<endl;
        } else {
        cout << "error :qpOASES solver failed due to infeasibility or other internal "
                "reasons" <<endl;
        }
    }

    //cout<<" reslut:"<<endl;
    double result[basic_control_size_*Nc_+1];
    qp_problem.getPrimalSolution(result);

    qp_problem.printOptions();
    optim_control_.resize(basic_control_size_*Nc_);
    //cout<<"size "<<basic_control_size_*Nc_<<endl;
    for(int i=0;i<basic_control_size_*Nc_;i++)
    {
        optim_control_[i]=result[i];
        //cout<<optim_control_[i]<<endl;
    }
        

    optim_relax_factor_=result[basic_control_size_*Nc_];

    //保存实际控制量，得到的最优控制增量+上一时刻的实际控制量
    real_control_vel_.push_back(optim_control_[0]+real_control_vel_.back());
    real_control_delta_.push_back(optim_control_[1]+real_control_delta_.back());
    
    return qp_problem.isSolved();
}


bool MPCControl::getFirstControl(double& vel,double& delta)
{
    if(optim_control_.size()<2)
    {
        cout<<"Error : not enough optim_control "<<endl;
        return false;
    }


    vel=real_control_vel_.back();
    delta=real_control_delta_.back();

    return true;
}


bool MPCControl::getRealControl(vector<double>& real_vel,vector<double>& real_delta)
{
    if(real_control_vel_.size()==0 | real_control_delta_.size()==0)
    {
        cout<<"MPCControl::getRealControl(vector<double>& real_vel,vector<double>& real_delta) :real_control size is zero "<<endl;
        return false;
    }

    real_vel.clear();
    real_delta.clear();
    real_vel.resize(real_control_vel_.size());
    real_delta.resize(real_control_delta_.size());
    copy(real_control_vel_.begin(),real_control_vel_.end(),real_vel.begin());
    copy(real_control_delta_.begin(),real_control_delta_.end(),real_delta.begin());

    for(int i=0;i<real_control_vel_.size();i++)
    {
      ; // cout<<real_control_vel_[i]<<" "<<real_vel[i]<<endl;

    }

    return true;
}


void MPCControl::updateErrorsRefTraj(vector<traj> & errors_ref_traj)
{
    errors_ref_traj_.resize(errors_ref_traj.size());
    for(int i=0;i<errors_ref_traj.size();i++)
    {
        errors_ref_traj_[i].x=errors_ref_traj[i].x-errors_ref_traj.front().x;
        errors_ref_traj_[i].y=errors_ref_traj[i].y-errors_ref_traj.front().y;
        errors_ref_traj_[i].phi=errors_ref_traj[i].phi-errors_ref_traj.front().phi;
        errors_ref_traj_[i].vel=errors_ref_traj[i].vel-errors_ref_traj.front().vel;
        errors_ref_traj_[i].delta=errors_ref_traj[i].delta-errors_ref_traj.front().delta;
    }
        
    
}