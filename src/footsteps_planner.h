#pragma once

#include <mc_control/mc_controller.h>

#include <ctime>

#include "polynomials.h"

namespace mc_plugin::footsteps_planner
{

struct Admissible_Region
{

public:
  Admissible_Region() = default;

  Admissible_Region(const sva::PTransformd & center, const Eigen::Vector2d & size)
  {
    compute_region(center.translation().segment(0, 2), mc_rbdyn::rpyFromMat(center.rotation()).z(), size);
  }
  Admissible_Region(const Eigen::Vector3d & center, const Eigen::Vector3d & size)
  {
    compute_region(center.segment(0, 2), center.z(), size.segment(0, 2));
  }
  void compute_region(const Eigen::Vector2d & center, double angle, const Eigen::Vector2d & size)
  {

    _center.segment(0, 2) = center;
    _angle = angle;
    _size.segment(0, 2) = size;
    _size.z() = 0;
    _center.z() = 0;
    R.setZero();
    R(0, 0) = cos(_angle);
    R(0, 1) = -sin(_angle);
    R(1, 0) = sin(_angle);
    R(1, 1) = cos(_angle);
    R(2, 2) = 1;
    upper_left_corner = _center + R * Eigen::Vector3d{-_size.x() / 2, _size.y() / 2, 0};
    upper_right_corner = _center + R * Eigen::Vector3d{_size.x() / 2, _size.y() / 2, 0};
    lower_left_corner = _center + R * Eigen::Vector3d{-_size.x() / 2, -_size.y() / 2, 0};
    lower_right_corner = _center + R * Eigen::Vector3d{_size.x() / 2, -_size.y() / 2, 0};

    std::vector<Eigen::Vector3d> Polygone_Corners = {upper_left_corner, upper_right_corner, lower_right_corner,
                                                     lower_left_corner};
    Polygone_Normals.resize(2, static_cast<Eigen::Index>(Polygone_Corners.size()));
    Polygone_Edges_Center.resize(2, Polygone_Normals.cols());
    Polygone_Vertices.resize(2, Polygone_Normals.cols());
    Offset.resize(Polygone_Normals.cols());
    for(size_t c = 0; c < Polygone_Corners.size(); c++)
    {
      Eigen::Vector3d point_1 = Polygone_Corners[c];
      Eigen::Vector3d point_2 = Polygone_Corners[(c + 1) % Polygone_Corners.size()];
      Eigen::Vector3d vertice = (point_2 - point_1).normalized();
      Eigen::Vector3d normal = Eigen::Vector3d{0, 0, 1}.cross(vertice).normalized();
      Polygone_Normals(0, static_cast<Eigen::Index>(c)) = normal.x();
      Polygone_Normals(1, static_cast<Eigen::Index>(c)) = normal.y();
      Polygone_Vertices(0, static_cast<Eigen::Index>(c)) = vertice.x();
      Polygone_Vertices(1, static_cast<Eigen::Index>(c)) = vertice.y();

      Eigen::Matrix2d R_Vertices_0;
      R_Vertices_0(0, 0) = Polygone_Normals(0, static_cast<Eigen::Index>(c));
      R_Vertices_0(1, 0) = Polygone_Vertices(0, static_cast<Eigen::Index>(c));
      R_Vertices_0(1, 0) = Polygone_Normals(1, static_cast<Eigen::Index>(c));
      R_Vertices_0(1, 1) = Polygone_Vertices(1, static_cast<Eigen::Index>(c));

      Offset(static_cast<Eigen::Index>(c)) = (R_Vertices_0.transpose() * Eigen::Vector2d{point_1.x(), point_1.y()}).x();
    }
  }
  ~Admissible_Region() = default;

  std::vector<Eigen::Vector3d> Get_corners()
  {
    return {upper_left_corner, upper_right_corner, lower_right_corner, lower_left_corner};
  }

  Eigen::MatrixXd Polygone_Normals;
  Eigen::VectorXd Offset;

private:
  Eigen::Vector3d _center = Eigen::Vector3d::Zero();
  Eigen::Vector3d _size = Eigen::Vector3d::Zero();
  double _angle = 0.0;
  Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
  Eigen::Vector3d upper_left_corner = Eigen::Vector3d::Zero();
  Eigen::Vector3d upper_right_corner = Eigen::Vector3d::Zero();
  Eigen::Vector3d lower_left_corner = Eigen::Vector3d::Zero();
  Eigen::Vector3d lower_right_corner = Eigen::Vector3d::Zero();
  Eigen::MatrixXd Polygone_Vertices;
  Eigen::MatrixXd Polygone_Edges_Center;
};

struct Steps_timings_output
{

public:
  Steps_timings_output() = default;
  ~Steps_timings_output() = default;
  bool QPsuccess;
  Eigen::VectorXd Ts;
  double loss;
};

struct ref_traj_point
{

  ref_traj_point() = default;
  ref_traj_point(const Eigen::Vector2d & pose, double ori)
  {
    pose_ = pose;
    ori_ = ori;
  }
  ref_traj_point(const sva::PTransformd & pose)
  {
    pose_ = pose.translation().segment(0, 2);
    ori_ = mc_rbdyn::rpyFromMat(pose.rotation()).z();
  }
  ~ref_traj_point() = default;

  const Eigen::Vector2d & pose() const noexcept
  {
    return pose_;
  }
  Eigen::Vector3d vec3_pose()
  {
    return Eigen::Vector3d{pose_.x(), pose_.y(), 0.};
  }
  sva::PTransformd PT_pose()
  {
    Eigen::Vector3d center{pose_.x(), pose_.y(), 0.};
    return sva::PTransformd(sva::RotZ(ori_), center);
  }
  double ori() const noexcept
  {
    return ori_;
  }
  Eigen::Vector3d rpy_ori() const noexcept
  {
    return Eigen::Vector3d{0., 0., ori_};
  }
  void ori(const double theta)
  {
    ori_ = theta;
  }

  Eigen::Vector2d pose_ = Eigen::Vector2d::Zero();
  double ori_ = 0.0;
};

struct Footstep : public mc_plugin::footsteps_planner::ref_traj_point
{

public:
  Footstep() = default;
  Footstep(const sva::PTransformd & pose, double ts, const Eigen::Vector2d & step_size)
  : mc_plugin::footsteps_planner::ref_traj_point(pose)
  {
    step_size_ = step_size;
    ts_ = ts;
    init();
  }
  Footstep(const Eigen::Vector2d & pose, const double ori, double ts, const Eigen::Vector2d & step_size)
  : mc_plugin::footsteps_planner::ref_traj_point(pose, ori)
  {
    step_size_ = step_size;
    ts_ = ts;
    init();
  }
  ~Footstep() = default;

  Admissible_Region & rect() noexcept
  {
    return step_rect_;
  }

  double ts() const noexcept
  {
    return ts_;
  }
  Eigen::Vector3d pose3()
  {
    return Eigen::Vector3d{pose_.x(), pose_.y(), 0};
  }
  const Eigen::Vector2d & size()
  {
    return step_size_;
  }

private:
  Eigen::Vector2d step_size_;
  Admissible_Region step_rect_;
  double ts_;

  inline void init() noexcept
  {
    Eigen::Vector3d center{pose_.x(), pose_.y(), ori_};
    Eigen::Vector3d dim{step_size_.x(), step_size_.y(), 0.};
    step_rect_ = Admissible_Region(center, dim);
  }
};

struct Footsteps_plan
{

public:
  Footsteps_plan() = default;
  Footsteps_plan(const Footstep & support_foot,
                 std::string_view supportFoot_name,
                 const Footstep & initial_swing_foot,
                 const std::vector<Footstep> & footsteps)
  {
    footsteps_ = footsteps;
    support_foot_name(supportFoot_name);
    support_foot_ = support_foot;
    initial_swing_foot_ = initial_swing_foot;
  }
  ~Footsteps_plan() = default;

  void add(Footstep step)
  {
    size_t i = 0;
    if(support_foot_name_ == "LeftFoot")
    {
      i = 1;
    }
    step.ori(step.ori() + ori_offset_ * std::pow(-1, footsteps_.size() + i));
    footsteps_.push_back(step);
  }
  void edit(const Footstep & step, const size_t indx)
  {
    footsteps_[indx] = step;
  }
  void clear()
  {
    footsteps_.clear();
  }
  std::vector<std::vector<Eigen::Vector3d>> get_steps_corners()
  {
    std::vector<std::vector<Eigen::Vector3d>> Output;
    for(size_t k = 0; k < n_steps(); k++)
    {
      Output.push_back(footsteps_[k].rect().Get_corners());
    }
    return Output;
  }
  const std::vector<Footstep> & steps()
  {
    return footsteps_;
  }
  std::vector<Eigen::Vector3d> steps_pose()
  {
    std::vector<Eigen::Vector3d> output;
    for(size_t k = 0; k < footsteps_.size(); k++)
    {

      output.push_back(footsteps_[k].vec3_pose());
    }
    return output;
  }
  std::vector<sva::PTransformd> steps_PTpose()
  {
    std::vector<sva::PTransformd> output;
    for(size_t k = 0; k < footsteps_.size(); k++)
    {

      output.push_back(footsteps_[k].PT_pose());
    }
    return output;
  }
  std::vector<double> steps_timings()
  {
    std::vector<double> output;
    for(size_t k = 0; k < footsteps_.size(); k++)
    {

      output.push_back(footsteps_[k].ts());
    }
    return output;
  }
  size_t n_steps() const noexcept
  {
    return footsteps_.size();
  }
  Footstep & support_foot()
  {
    return support_foot_;
  }
  Footstep & footstep(size_t indx)
  {
    return footsteps_[indx];
  }
  void support_foot(const Footstep & footstep)
  {
    support_foot_no_offset = footstep;
    support_foot_ = footstep;
    if((footstep.pose() - support_foot_.pose()).norm() > 1e-3 && std::abs(footstep.ori() - support_foot_.ori()) > 1e-3
       && !offset_applied)
    {
      offset_applied = true;
    }
    if(offset_applied)
    {
      double sgn = -1;
      if(support_foot_name_ == "LeftFoot")
      {
        sgn *= -1;
      }
      support_foot_.ori(footstep.ori() - sgn * ori_offset_);
    }
  }

  double ori_offset() const noexcept
  {
    return ori_offset_;
  }

  void ori_offset(const double & theta)
  {
    if(theta != ori_offset_)
    {
      offset_applied = false;
    }
    ori_offset_ = theta;
  }

  void support_foot_name(std::string_view name) noexcept
  {
    support_foot_name_ = name;
  }

  const std::string & support_foot_name()
  {
    return support_foot_name_;
  }

private:
  std::vector<Footstep> footsteps_;
  Footstep support_foot_;
  Footstep support_foot_no_offset;
  std::string support_foot_name_ = "LeftFoot";
  Footstep initial_swing_foot_;
  double ori_offset_ = 0;
  bool offset_applied = false;
};

struct FootStepGen
{

public:
  FootStepGen(const mc_rtc::Configuration & config);
  FootStepGen();

  ~FootStepGen() = default;
  /**
   * Initialize the footsteps Generator
   * @tparam supportFootName
   * @tparam P_f0 Support Foot
   * @tparam V Reference velocity inputs from t0 (can be empty)
   * @tparam Tstep desired ordered Steps Timing (can be empty)
   * @tparam Pf desired ordered Footsteps coordinate (can be empty) (angle in z coordinate)
   */
  void init(std::string_view supportFootName,
            const Footstep & P_f0,
            const std::vector<sva::MotionVecd> & V,
            const std::vector<double> & Tstep,
            const std::vector<Footstep> & ref_pose);

  void reconfigure(const mc_rtc::Configuration & config);

  // Return The footsteps Theta values
  const Eigen::VectorXd & Theta_f() const noexcept
  {
    return Theta_f_;
  }

  // Return the reference trajectory in the preview horizon
  std::vector<Eigen::Vector3d> ref_traj(bool centered)
  {
    std::vector<Eigen::Vector3d> Output;

    for(size_t k = 0; k < P_traj_.size(); k++)
    {
      Eigen::Vector3d offset = P_traj_[0].vec3_pose();
      if(centered)
      {
        offset = Eigen::Vector3d::Zero();
      }
      Output.push_back(P_traj_[k].vec3_pose() - offset);
    }
    return Output;
  }
  // Compute The Footsteps and the Steps Timings
  Footsteps_plan compute_plan();

  Footsteps_plan & footsteps_plan() noexcept
  {
    return plan_;
  }

  Eigen::Index Get_Nsteps() const noexcept
  {
    return N_steps;
  }

  Eigen::Vector2d intersec = Eigen::Vector2d::Zero();
  Eigen::Vector2d r = Eigen::Vector2d::Zero();

private:
  Eigen::VectorXd solveQP();

  /**
   * Compute N points trajectory between P_s_0 and P_s_1
   * @return Points Coordonate and angle of the trajectory
   */
  std::vector<ref_traj_point> GetRefTrajectory(ref_traj_point & P_s_0, ref_traj_point & P_s_1, double duration);

  std::vector<std::vector<double>> GetVelocityProfile(const Eigen::Vector3d & P_s_0,
                                                      double V_Max,
                                                      double V_Min,
                                                      const std::vector<Eigen::Vector3d> & Traj);

  // Compute the Steps Timing dependings of the given parameter
  void GetStepsTimings();

  /**
   * return the position of the reference velocity integratin the velocity profile
   * @tparam k_end time index desired
   * @return Coordinate of the integrated ref velocity at time index k with orientation in z
   */
  ref_traj_point IntegrateVelProfile(size_t k_end);

  std::string supportFoot = "RightFoot";

  HoubaPolynomial<Eigen::Vector2d> path;

  Footsteps_plan plan_;
  std::vector<Footstep> pose_reference_;
  std::vector<Footstep> steps_inputs_;
  std::vector<double> t_steps_inputs_; // Input Step Timing
  std::vector<sva::MotionVecd> v_inputs_; // Velocity input

  Eigen::Index F_ = 1; // footsteps number
  Eigen::Index N_steps = -1;

  std::vector<double> StepsTimings_; // Contains the time of each steps
  std::vector<size_t> StepsTimings_indx_; // Index timing of each steps
  std::vector<int> FootSteps_indx_; // Index of the input steps position for the right step timing
  Eigen::VectorXd Theta_f_; // Output Steps Angles
  Eigen::VectorXd m_Ts; // Steps Duration

  std::vector<ref_traj_point> P_traj_; // Position of reference trajectory for each timesteps

  // QP Problem
  bool QPsuccess = false;
  Eigen::MatrixXd Q_; // QP Hessian
  Eigen::VectorXd p_; // QP Grad

  Eigen::MatrixXd Aeq; // Equality Matrix
  Eigen::VectorXd beq; // Equality Vector

  Eigen::MatrixXd Aineq; // Inequality Matrix
  Eigen::VectorXd bineq; // Inequality Vector

  bool legs_warning_on = false;

public:
  double Ts_min_ = 0.8; // Step Time lenght limits
  double Ts_max_ = 2; // Step Time lenght limits
  double l_ = 0.2; // Distance between foot
  double d_min = 0.34; // min distance between the foot center and the kinematic constraint;
  double Tp_ = 6; // Preview horizon time
  double delta_ = 5e-2; // t_k - t_k-1
  double d_h_x = 0.2; // Next step tolerance zone
  double d_h_y = 0.05; // Next step tolerance zone
  double v_ = 0.1; // Cruise Parameters
  double max_theta = 3.14 / 6; // Max angle between two steps
  int P_ = 100; // Preview horizon time indexes
  double Ts_ = 5.0; // Cruise Parameters
  double robot_height_ = 150; // in cm
  double theta_offset_ = 0;

  Eigen::Vector3d v_ref_ = Eigen::Vector3d::Zero();
};

} // namespace mc_plugin::footsteps_planner
