#include "footsteps_planner.h"

#include <eigen-quadprog/QuadProg.h>

namespace mc_plugin::footsteps_planner
{

void FootStepGen::init(std::string_view supportFootName,
                       const Footstep & P_f0,
                       const std::vector<sva::MotionVecd> & V,
                       const std::vector<double> & Tstep,
                       const std::vector<Footstep> & ref_pose)
{
  plan_.support_foot_name(supportFootName);
  plan_.support_foot(P_f0);
  N_steps = -1;
  supportFoot = supportFootName;
  v_inputs_ = V;
  P_traj_.clear();
  P_ = static_cast<int>(Tp_ / delta_);
  for(size_t k = 0; k < v_inputs_.size(); k++)
  {
    P_traj_.push_back(IntegrateVelProfile(k));
  }
  t_steps_inputs_ = Tstep;

  if(!t_steps_inputs_.empty())
  {
    t_steps_inputs_[0] = std::max(Ts_min_, t_steps_inputs_[0]);
  }

  pose_reference_ = ref_pose;
  for(auto & ref : pose_reference_)
  {
    if(std::abs(ref.ori()) > M_PI)
    {
      ref.ori(ref.ori() - (ref.ori() / std::abs(ref.ori())) * 2 * M_PI);
    }
  }
}

void FootStepGen::GetStepsTimings()
{
  StepsTimings_.clear();
  StepsTimings_ = t_steps_inputs_;
  StepsTimings_indx_.clear();
  FootSteps_indx_.clear();
  // 1 Generate reference trajectory and velocity according to order V -> Pf -> Ts

  if(!pose_reference_.empty())
  {

    ref_traj_point P_f_im1 = IntegrateVelProfile(0);
    ref_traj_point P_f_i(pose_reference_[0].pose(), pose_reference_[0].ori());
    const double d_step = Eigen::Vector2d{d_h_x, d_h_y}.norm() / 2;
    const size_t nsteps = static_cast<size_t>((P_f_im1.pose() - P_f_i.pose()).norm() / d_step);
    if(nsteps < StepsTimings_.size())
    {
      P_traj_ = GetRefTrajectory(P_f_im1, P_f_i, StepsTimings_[nsteps]);
      while(P_traj_.size() < static_cast<size_t>(P_))
      {
        P_traj_.push_back(P_traj_.back());
      }
    }
    else
    {
      P_traj_ = GetRefTrajectory(P_f_im1, P_f_i, Tp_);
    }
  }

  for(size_t i = 0; i < StepsTimings_.size(); i++)
  {
    double t_i = StepsTimings_[i];
    StepsTimings_indx_.push_back(static_cast<size_t>(std::round(t_i / delta_)));
    FootSteps_indx_.push_back(-1);
  }
  while(StepsTimings_[StepsTimings_.size() - 1] > Tp_ - delta_)
  {
    StepsTimings_.pop_back();
    FootSteps_indx_.pop_back();
    StepsTimings_indx_.pop_back();
  }

  F_ = static_cast<int>(StepsTimings_.size());
}

Footsteps_plan FootStepGen::compute_plan()
{
  plan_.clear();
  GetStepsTimings();

  Theta_f_.resize(F_, 1);

  if(F_ <= 0)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>("[footstepsplanner::compute_plan] No step to compute");
    return plan_;
  }

  // Solving QP 1 for Orientation

  Eigen::VectorXd Dtheta_upper_lim(Eigen::VectorXd::Ones(F_) * max_theta);
  Eigen::VectorXd Dtheta_lower_lim(Eigen::VectorXd::Ones(F_) * -max_theta);
  Eigen::MatrixXd Delta =
      Eigen::MatrixXd::Identity(F_, F_); // Differenciation matrix between two steps orientation (for constraints)
  Aeq = Eigen::MatrixXd::Zero(F_, F_);
  beq = Eigen::VectorXd::Zero(Aeq.rows());

  for(Eigen::Index k = 0; k < F_; ++k)
  {
    size_t k_sz = static_cast<size_t>(k);
    if(k == 0)
    {
      Dtheta_upper_lim(k) += plan_.support_foot().ori();
      Dtheta_lower_lim(k) += plan_.support_foot().ori();
    }

    else
    {

      Delta(k, k - 1) = -1;
    }

    if(FootSteps_indx_[k_sz] != -1)
    {
      size_t fs_idx = static_cast<size_t>(FootSteps_indx_[k_sz]);
      double theta_cstr = pose_reference_[fs_idx].ori();
      if(P_traj_[StepsTimings_indx_[k_sz]].ori() - pose_reference_[fs_idx].ori() > M_PI)
      {
        theta_cstr -= pose_reference_[fs_idx].ori() / std::abs(pose_reference_[fs_idx].ori()) * 2 * M_PI;
      }
      Aeq(k, k) = 1;
      beq(k) = theta_cstr;
    }
  }

  Aineq = Eigen::MatrixXd::Zero(2 * Delta.rows(), F_);
  bineq = Eigen::VectorXd::Zero(Aineq.rows());
  Aineq << Delta, -Delta;
  bineq << Dtheta_upper_lim, -Dtheta_lower_lim;

  Eigen::VectorXd b = Eigen::VectorXd::Zero(F_);
  for(Eigen::Index i = 0; i < F_; i++)
  {
    b(i) = P_traj_[StepsTimings_indx_[static_cast<size_t>(i)]].ori();
  }
  Q_ = Eigen::MatrixXd::Identity(F_, F_);
  p_ = -b;
  Eigen::VectorXd theta = solveQP();
  if(!QPsuccess)
  {
    mc_rtc::log::error("[Footsteps planner] Step theta failed");
  }

  for(int k = 0; k < theta.size(); k++)
  {

    if(std::abs(theta(k)) > M_PI)
    {
      Theta_f_(k) = theta(k) - (theta(k) / std::abs(theta(k))) * 2 * M_PI;
    }
    else
    {
      Theta_f_(k) = theta(k);
    }
  }

  // Solving QP 2 For placement

  Delta =
      Eigen::MatrixXd::Identity(2 * F_, 2 * F_); // Differenciation matrix between two steps location (for constraints)
  Aeq = Eigen::MatrixXd::Zero(2 * F_, 2 * F_);
  beq = Eigen::VectorXd::Zero(2 * F_);

  double sgn = -1.0;
  if(supportFoot == "RightFoot")
  {
    sgn = 1.0;
  }
  double l = l_;
  if(d_h_y / 2 > l_ - 0.5 * d_min)
  {
    l = 0.5 * (d_h_y + d_min);
  }

  std::vector<Eigen::VectorXd> cstr_vec;
  std::vector<Eigen::MatrixX2d> Normal_Vec;

  Eigen::Matrix2d R_Theta_0;
  double theta_0 = P_traj_[0].ori();
  R_Theta_0 = sva::RotZ(-theta_0).block(0, 0, 2, 2);

  Admissible_Region Kinematic_Rectangle(Eigen::Vector3d{0, 0, plan_.support_foot().ori()},
                                        Eigen::Vector3d{d_h_x, d_h_y, 0});

  Eigen::VectorXd bcstr = Kinematic_Rectangle.Offset
                          + Kinematic_Rectangle.Polygone_Normals.transpose()
                                * (plan_.support_foot().pose() + sgn * R_Theta_0 * Eigen::Vector2d{0, l});

  Normal_Vec.push_back(Kinematic_Rectangle.Polygone_Normals.transpose());
  cstr_vec.push_back(bcstr);

  for(Eigen::Index k = 1; k < F_; k++)
  {
    double theta_k = P_traj_[StepsTimings_indx_[static_cast<size_t>(k)]].ori();
    double theta_f_km1 = Theta_f_(k - 1);

    Eigen::Matrix2d R_Theta_k;
    R_Theta_k = sva::RotZ(-theta_k).block(0, 0, 2, 2);

    Kinematic_Rectangle = Admissible_Region(Eigen::Vector3d{0, 0, theta_f_km1}, Eigen::Vector3d{d_h_x, d_h_y, 0});
    Eigen::Matrix2d R_Theta_f_km1 = sva::RotZ(-theta_f_km1).block(0, 0, 2, 2);

    bcstr = Kinematic_Rectangle.Offset
            + (k % 2 == 0 ? sgn : -sgn) * Kinematic_Rectangle.Polygone_Normals.transpose() * R_Theta_f_km1
                  * Eigen::Vector2d{0, l};

    Normal_Vec.push_back(Kinematic_Rectangle.Polygone_Normals.transpose());
    cstr_vec.push_back(bcstr);
  }

  if(FootSteps_indx_[0] != -1)
  {
    size_t fs_idx = static_cast<size_t>(FootSteps_indx_[0]);
    Aeq(0, 0) = 1;
    Aeq(1, 1) = 1;
    beq(0) = steps_inputs_[fs_idx].pose().x();
    beq(1) = steps_inputs_[fs_idx].pose().y();
  }
  for(Eigen::Index k = 1; k < F_; k++)
  {

    Delta(2 * k, 2 * (k - 1)) = -1;
    Delta(2 * k + 1, 2 * (k - 1) + 1) = -1;

    if(FootSteps_indx_[static_cast<size_t>(k)] != -1)
    {
      size_t fs_idx = static_cast<size_t>(FootSteps_indx_[static_cast<size_t>(k)]);
      Aeq(2 * k, 2 * k) = 1;
      Aeq(2 * k + 1, 2 * k + 1) = 1;
      beq(2 * k) = steps_inputs_[fs_idx].pose().x();
      beq(2 * k + 1) = steps_inputs_[fs_idx].pose().y();
    }
  }

  Eigen::Index N_cstr = 0;
  for(size_t k = 0; k < static_cast<size_t>(Normal_Vec.size()); k++)
  {
    N_cstr += Normal_Vec[k].rows();
  }

  Aineq = Eigen::MatrixXd::Zero(N_cstr, 2 * F_);
  bineq = Eigen::VectorXd::Zero(N_cstr);

  Eigen::Index step_indx = 0;
  Eigen::Index cstr_index = 0;
  for(size_t i_ineq = 0; i_ineq < Normal_Vec.size(); i_ineq++)
  {

    Eigen::MatrixXd n_vec = Normal_Vec[i_ineq];
    Eigen::VectorXd ineq = cstr_vec[i_ineq];

    for(Eigen::Index cstr = 0; cstr < n_vec.rows(); cstr++)
    {

      Aineq(cstr_index + cstr, step_indx) = n_vec(cstr, 0);
      Aineq(cstr_index + cstr, step_indx + 1) = n_vec(cstr, 1);

      bineq(cstr_index + cstr) = ineq(cstr);
    }

    step_indx += 2;
    cstr_index += n_vec.rows();
  }

  Aineq = Aineq * Delta;

  b = Eigen::VectorXd::Zero(2 * F_);
  for(Eigen::Index i = 0; i < F_; i++)
  {
    size_t i_sz = static_cast<size_t>(i);
    double theta_i = P_traj_[i_sz].ori();
    Eigen::Matrix2d R_Theta_i;
    R_Theta_i = sva::RotZ(-theta_i).block(0, 0, 2, 2);
    Eigen::Vector2d dl = (i % 2 == 0 ? sgn : -sgn) * R_Theta_i * Eigen::Vector2d{0, l_ / 2};
    b(2 * i) = P_traj_[StepsTimings_indx_[i_sz]].pose().x() + dl.x();
    b(2 * i + 1) = P_traj_[StepsTimings_indx_[i_sz]].pose().y() + dl.y();
  }

  Q_ = Eigen::MatrixXd::Identity(2 * F_, 2 * F_);
  p_ = -b;
  Eigen::VectorXd XY(solveQP());
  if(!QPsuccess)
  {
    mc_rtc::log::error("[Footsteps planner] Step QP failed");
  }

  plan_.ori_offset(theta_offset_);

  for(Eigen::Index k = 0; k < F_; k++)
  {
    double xf = XY(2 * k);
    double yf = XY(2 * k + 1);
    plan_.add(Footstep(Eigen::Vector2d{xf, yf}, Theta_f_(k), StepsTimings_[static_cast<size_t>(k)],
                       Eigen::Vector2d{0.1, 0.1}));
  }

  const double theta_legs_0 = plan_.support_foot().ori();
  const double theta_legs_1 = Theta_f_(0);

  Eigen::Matrix2d A;
  A.block(0, 0, 2, 1) << cos(theta_legs_0), sin(theta_legs_0);
  A.block(0, 1, 2, 1) << -cos(theta_legs_1), -sin(theta_legs_1);
  if(A.determinant() != 0)
  {
    const Eigen::Matrix2d R_sup_0 = plan_.support_foot().PT_pose().rotation().block(0, 0, 2, 2).transpose();
    const Eigen::Vector2d coeff = A.inverse() * (plan_.steps_pose()[0].segment(0, 2) - plan_.support_foot().pose());
    intersec = plan_.support_foot().pose() + R_sup_0 * Eigen::Vector2d{1, 0} * coeff(0);
    const double proj = (intersec - plan_.support_foot().pose()).transpose()
                        * (plan_.steps_pose()[0].segment(0, 2) - plan_.support_foot().pose()).normalized();

    r = (intersec - plan_.support_foot().pose())
        - (plan_.steps_pose()[0].segment(0, 2) - plan_.support_foot().pose()).normalized() * proj;
    if((R_sup_0.transpose() * r).x() > 0 && r.norm() < 1)
    {
      if(!legs_warning_on)
      {
        mc_rtc::log::warning("[Footsteps planner] Potential Legs collision");
        legs_warning_on = true;
      }
      plan_.edit(
          Footstep(plan_.steps()[0].pose_, plan_.support_foot().ori(), StepsTimings_[0], Eigen::Vector2d{0.1, 0.1}), 0);
    }
    else
    {
      legs_warning_on = false;
    }
  }

  return plan_;
}

std::vector<ref_traj_point> FootStepGen::GetRefTrajectory(ref_traj_point & P_s_0,
                                                          ref_traj_point & P_s_1,
                                                          double duration)
{
  const Eigen::Matrix2d R_s0_0 = P_s_0.PT_pose().rotation().block(0, 0, 2, 2).transpose();

  bool shuffle = false;
  bool backward = false;

  if(std::abs((R_s0_0.transpose() * (P_s_1.pose() - P_s_0.pose())).x()) < 3e-1)
  {
    shuffle = true;
  }
  if((R_s0_0.transpose() * (P_s_1.pose() - P_s_0.pose())).x() < 0 && !shuffle)
  {
    backward = true;
  }

  const Eigen::Vector2d & init_pose = P_s_0.pose();
  const Eigen::Vector2d & target_pose = P_s_1.pose();

  if(!shuffle)
  {
    Eigen::Vector2d init_ori = {cos(P_s_0.ori()), sin(P_s_0.ori())};
    Eigen::Vector2d target_ori = {cos(P_s_1.ori()), sin(P_s_1.ori())};
    if(backward)
    {
      init_ori = {-cos(P_s_0.ori()), -sin(P_s_0.ori())};
      target_ori = {-cos(P_s_1.ori()), -sin(P_s_1.ori())};
    }
    path.reset(init_pose, init_ori, target_pose, target_ori);
  }

  const int N = static_cast<int>(duration / delta_);

  std::vector<ref_traj_point> Output;
  for(int k = 0; k < N + 1; k++)
  {
    double t = ((double)k) / ((double)N);
    if(!shuffle)
    {
      Eigen::Vector2d pos_t = path.pos(t);
      Eigen::Vector2d ori_t = path.tangent(t);
      double theta = atan2(ori_t.y(), ori_t.x());
      if((init_pose - target_pose).norm() < 5e-2)
      {
        theta = P_s_1.ori();
      }
      if(backward)
      {
        const double sgn = theta / std::abs(theta);
        Output.push_back(ref_traj_point(pos_t, std::fmod(theta - sgn * M_PI, 2 * M_PI)));
      }
      else
      {
        Output.push_back(ref_traj_point(pos_t, std::fmod(theta, 2 * M_PI)));
      }
    }
    else
    {
      Output.push_back(ref_traj_point(P_s_0.pose() + t * (P_s_1.pose() - P_s_0.pose()), P_s_1.ori()));
    }
  }

  return Output;
}

Eigen::VectorXd FootStepGen::solveQP()
{

  Eigen::QuadProgDense QP;

  int Nvar = static_cast<int>(Q_.rows());
  int NIneqConstr = static_cast<int>(Aineq.rows());
  int NEqConstr = static_cast<int>(Aeq.rows());
  QP.problem(Nvar, NEqConstr, NIneqConstr);
  QPsuccess = QP.solve(Q_, p_, Aeq, beq, Aineq, bineq);

  return QP.result();
}

ref_traj_point FootStepGen::IntegrateVelProfile(size_t k_end)
{

  size_t k_start = 0;
  if(k_end < P_traj_.size())
  {
    return P_traj_[k_end];
  }

  ref_traj_point Output(plan_.support_foot().pose(), plan_.support_foot().ori());
  if(P_traj_.empty())
  {

    if(supportFoot == "RightFoot")
    {
      Output.pose_ += sva::RotZ(-Output.ori_).block(0, 0, 2, 2) * Eigen::Vector2d{0, l_ / 2};
    }
    else
    {
      Output.pose_ -= sva::RotZ(-Output.ori_).block(0, 0, 2, 2) * Eigen::Vector2d{0, l_ / 2};
    }
  }
  else
  {
    Output.pose_ = P_traj_.back().pose();
    Output.ori_ = P_traj_.back().ori();
    k_start = P_traj_.size() - 1;
  }

  for(size_t i = k_start; i < std::min(k_end, v_inputs_.size()); i++)
  {
    Output.ori_ += v_inputs_[i].angular().z() * delta_;
    Output.pose_ += sva::RotZ(-Output.ori_).block(0, 0, 2, 2) * v_inputs_[i].linear().segment(0, 2) * delta_;
  }

  return Output;
}

} // namespace mc_plugin::footsteps_planner
