#ifndef MATRIX_H
#define MATRIX_H
#include "eigen3/Eigen/Dense"

namespace Oceane {

using Matrix=  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> ;
using Vector= Eigen::VectorXd;

} // namespace Oceane

#endif // MATRIX_H
