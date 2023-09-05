#include "quad4.h"

namespace Oceane {


/*.........................................................
 Implementation of Global class methods..i.e CPS4..........
...........................................................*/

void Global::getBsmat_at(double x, double y, Oceane::Matrix &Bs, double &detJ)
{
    if(Bs.rows() !=3 or Bs.cols() != 24 ) Bs.resize(3,24);
    Eigen::Matrix<double,24,24> coeff;
    for(long i=0; i<m_xvec.size(); ++i)
    {

        auto  xval= m_xvec[i]; auto yval =m_yvec[i];
        Eigen::Matrix<double,6,24> temp;
        temp<<
               1 ,xval, yval, pow(xval,2), pow(yval,2), xval*yval, pow(xval,3), pow(xval,2)*yval, xval*pow(yval,2), pow(yval,3), pow(xval,4), pow(xval,3)*yval, pow(xval,2)*pow(yval,2), xval*pow(yval,3), pow(yval,4), pow(xval,4)*yval, pow(xval,3)*pow(yval,2), pow(xval,2)*pow(yval,3), xval*pow(yval,4), pow(xval,5)*yval, pow(xval,3)*pow(yval,3), xval*pow(yval,5), pow(xval,5), pow(yval,5),
               0, 1, 0, 2*xval, 0, yval, 3*pow(xval,2), 2*xval*yval, pow(yval,2), 0, 4*pow(xval,3), 3*pow(xval,2)*yval, 2*xval*pow(yval,2), pow(yval,3), 0, 4*pow(xval,3)*yval, 3*pow(xval,2)*pow(yval,2), 2*xval*pow(yval,3), pow(yval,4), 5*pow(xval,4)*yval, 3*pow(xval,2)*pow(yval,3), pow(yval,5), 5*pow(xval,4), 0,
               0, 0, 1, 0, 2*yval, xval, 0, pow(xval,2), 2*xval*yval, 3*pow(yval,2), 0, pow(xval,3), 2*pow(xval,2)*yval, 3*xval*pow(yval,2), 4*pow(yval,3), pow(xval,4), 2*pow(xval,3)*yval, 3*pow(xval,2)*pow(yval,2), 4*xval*pow(yval,3), pow(xval,5), 3*pow(xval,3)*pow(yval,2), 5*xval*pow(yval,4), 0, 5*pow(yval,4),
               0, 0, 0, 2, 0, 0, 6*xval, 2*yval, 0, 0, 12*pow(xval,2), 6*xval*yval, 2*pow(yval,2), 0, 0, 12*pow(xval,2)*yval, 6*xval*pow(yval,2), 2*pow(yval,3), 0, 20*pow(xval,3)*yval, 6*xval*pow(yval,3), 0, 20*pow(xval,3), 0,
               0, 0, 0, 0, 2, 0, 0, 0, 2*xval, 6*yval, 0, 0, 2*pow(xval,2), 6*xval*yval, 12*pow(yval,2), 0, 2*pow(xval,3), 6*pow(xval,2)*yval, 12*xval*pow(yval,2), 0, 6*pow(xval,3)*yval, 20*xval*pow(yval,3), 0, 20*pow(yval,3),
               0, 0, 0, 0, 0, 1, 0, 2*xval, 2*yval, 0, 0, 3*pow(xval,2), 4*xval*yval, 3*pow(yval,2), 0, 4*pow(xval,3), 6*pow(xval,2)*yval, 6*xval*pow(yval,2), 4*pow(yval,3), 5*pow(xval,4), 9*pow(xval,2)*pow(yval,2), 5*pow(yval,4), 0, 0;

        coeff.block<6,24>(i*6,0) = temp;
    }

    Eigen::Matrix<double,24,24> I;
    I.setIdentity();
    Eigen::MatrixXd A = coeff.fullPivLu().solve(I);
    A.transposeInPlace();

    Eigen::Matrix<double,24,1> phi_yy_coeff;
    phi_yy_coeff<<  0, 0, 0, 0, 2, 0, 0, 0, 2*x, 6*y, 0, 0, 2*pow(x,2), 6*x*y, 12*pow(y,2), 0, 2*pow(x,3), 6*pow(x,2)*y, 12*x*pow(y,2), 0, 6*pow(x,3)*y, 20*x*pow(y,3), 0, 20*pow(y,3);
    Eigen::Matrix<double,24,1> phi_xx_coeff;
    phi_xx_coeff<<  0, 0, 0, 2, 0, 0, 6*x, 2*y, 0, 0, 12*pow(x,2), 6*x*y, 2*pow(y,2), 0, 0, 12*pow(x,2)*y, 6*x*pow(y,2), 2*pow(y,3), 0, 20*pow(x,3)*y, 6*x*pow(y,3), 0, 20*pow(x,3), 0;
    Eigen::Matrix<double,24,1> phi_xy_coeff;
    phi_xy_coeff<<  0, 0, 0, 0, 0, -1, 0, -2*x, -2*y, 0, 0, -3*pow(x,2), -4*x*y, -3*pow(y,2), 0, -4*pow(x,3), -6*pow(x,2)*y, -6*x*pow(y,2), -4*pow(y,3), -5*pow(x,4), -9*pow(x,2)*pow(y,2), -5*pow(y,4), 0, 0;


    assert(Bs.rows()==3 and Bs.cols()==24);
    Bs.row(0) = (A*phi_yy_coeff).transpose();
    Bs.row(1) = (A*phi_xx_coeff).transpose();
    Bs.row(2) = (A*phi_xy_coeff).transpose();

    detJ =1.0;
}


void Global::getBdmat_at(double x, double y, Oceane::Matrix& bd, double& detJ)
{

    if(bd.rows()!=3 or bd.cols()!=8)
        bd.resize(3,8);
    Eigen::Matrix<double,4,4> coeff;
    coeff<<
            1,m_xvec[0],m_yvec[0],m_xvec[0]*m_yvec[0],
            1,m_xvec[1],m_yvec[1],m_xvec[1]*m_yvec[1],
            1,m_xvec[2],m_yvec[2],m_xvec[2]*m_yvec[2],
            1,m_xvec[3],m_yvec[3],m_xvec[3]*m_yvec[3] ;
    assert(coeff.determinant()!=0.0);
    Oceane::Matrix disp_coeff = coeff.inverse().transpose();


   Oceane::Vector dudx(4),dudy(4);
    dudx[0]=0; dudx[1]=1; dudx[2]=0; dudx[3]=y;
    dudy[0]=0; dudy[1]=0; dudy[2]=1; dudy[3]=x;

    bd.setZero(); Oceane::Vector dx,dy;
    dx = disp_coeff*dudx;
    dy = disp_coeff*dudy;

    bd(0,0)=dx[0]; bd(0,2)=dx[1]; bd(0,4)=dx[2]; bd(0,6)=dx[3];
    bd(1,1)=dy[0]; bd(1,3)=dy[1]; bd(1,5)=dy[3]; bd(1,7)=dy[3];
    bd(2,0)=dy[0]; bd(2,2)=dy[1]; bd(2,4)=dy[2]; bd(2,6)=dy[3];
    bd(2,1)=dx[0]; bd(2,3)=dx[1]; bd(2,5)=dx[2]; bd(2,7)=dx[3];

    detJ =1.0;
}









/*.........................................................
..................PARAMETRIC CLASS DEFINITIONS.............
 *........................................................*/
void Shape::addCoordinate(Oceane::Vector xcoord, Oceane::Vector ycoord)
{
    m_xcoord=xcoord;
    m_ycoord=ycoord;
}

void Shape::update(double sval, double tval)
{
    m_dNx[0]=(tval-1)/4;
    m_dNx[1]=-(tval-1)/4;
    m_dNx[2]=(tval+1)/4;
    m_dNx[3]=-(tval+1)/4;

    m_dNy[0]= (sval-1)/4;
    m_dNy[1]=-(sval+1)/4;
    m_dNy[2]=(sval+1)/4;
    m_dNy[3]=-(sval-1)/4;

    m_JacobMatrix(0,0) = m_dNx.transpose() * m_xcoord;
    m_JacobMatrix(0,1) = m_dNy.transpose() * m_xcoord;
    m_JacobMatrix(1,0) = m_dNx.transpose() * m_ycoord;
    m_JacobMatrix(1,1) = m_dNy.transpose() * m_ycoord;

}

Oceane::Matrix Shape::getJacobMatrix()
{
    return m_JacobMatrix;
}

Oceane::Vector Shape::getdNx()
{
    return m_dNx;
}

Oceane::Vector Shape::getdNy()
{
    return m_dNy;
}



void Parametric::addCoordinates(Oceane::Vector xvec,Oceane::Vector yvec)
{
    m_shape.addCoordinate(xvec,yvec);
}



void Parametric::getBdmat_at(double x, double y, Oceane::Matrix& Bdmat, double& detJ)
{
    if(Bdmat.rows()!=3 or (Bdmat.cols()!=8)) Bdmat.resize(3,8);

    m_shape.update(x,y);

    detJ= m_shape.getJacobMatrix().determinant();
    Eigen::MatrixXd inv= m_shape.getJacobMatrix().inverse();
    Eigen::VectorXd Nis,Nit;
    Nis=m_shape.getdNx();
    Nit=m_shape.getdNy();
    Bdmat.setZero();

    Bdmat(0,0)=Nis(0)*inv(0,0)+Nit(0)*inv(1,0);
    Bdmat(0,2)=Nis(1)*inv(0,0)+Nit(1)*inv(1,0);
    Bdmat(0,4)=Nis(2)*inv(0,0)+Nit(2)*inv(1,0);
    Bdmat(0,6)=Nis(3)*inv(0,0)+Nit(3)*inv(1,0);

    Bdmat(1,1)=Nis(0)*inv(0,1)+Nit(0)*inv(1,1);
    Bdmat(1,3)=Nis(1)*inv(0,1)+Nit(1)*inv(1,1);
    Bdmat(1,5)=Nis(2)*inv(0,1)+Nit(2)*inv(1,1);
    Bdmat(1,7)=Nis(3)*inv(0,1)+Nit(3)*inv(1,1);

    Bdmat(2,0)=Nis(0)*inv(0,1)+Nit(0)*inv(1,1);
    Bdmat(2,2)=Nis(1)*inv(0,1)+Nit(1)*inv(1,1);
    Bdmat(2,4)=Nis(2)*inv(0,1)+Nit(2)*inv(1,1);
    Bdmat(2,6)=Nis(3)*inv(0,1)+Nit(3)*inv(1,1);

    Bdmat(2,1)=Nis(0)*inv(0,0)+Nit(0)*inv(1,0);
    Bdmat(2,3)=Nis(1)*inv(0,0)+Nit(1)*inv(1,0);
    Bdmat(2,5)=Nis(2)*inv(0,0)+Nit(2)*inv(1,0);
    Bdmat(2,7)=Nis(3)*inv(0,0)+Nit(3)*inv(1,0);

}


void Parametric::getBsmat_at(double x, double y, Oceane::Matrix& Bs, double& detJ)
{

    m_shape.update(x,y);
    detJ = m_shape.getJacobMatrix().determinant();

    Bs.setZero();
    Bs(0,0)=(3*y*(x - 1)*(pow(x,2) + x + 5*pow(y,2) - 5))/8;
    Bs(0,1)= 3*y*pow(x -1 , 2)*(x + 1)/8;
    Bs(0,2) = (x - 1)*(3*pow(x , 2)*y - pow(x , 2) + 3*x*y - x + 15*pow(y,3) - 3*pow(y,2) - 15*y + 3)/8;
    Bs(0,4) = (x - 1)*(y - 1)*(5*pow(y,2) + 2*y - 1)/8;
    Bs(0,5) = (3*y - 1)*pow(x -1 , 2)*(x + 1)/8;
    Bs(0,6) = 3*y*(x + 1)*(- pow(x , 2) + x - 5*pow(y,2) + 5)/8;
    Bs(0,7) = 3*y*(x - 1)*pow(x +1 , 2)/8;
    Bs(0,8) = -(x + 1)*(3*pow(x , 2)*y - pow(x , 2) - 3*x*y + x + 15*pow(y,3) - 3*pow(y,2) - 15*y + 3)/8;
    Bs(0,10) = -(x + 1)*(y - 1)*(5*pow(y,2) + 2*y - 1)/8;
    Bs(0,11) = (3*y - 1)*(x - 1)*pow(x +1 , 2)/8;
    Bs(0,12) = 3*y*(x + 1)*(pow(x , 2) - x + 5*pow(y,2) - 5)/8;
    Bs(0,13) = -3*y*(x - 1)*pow(x +1 , 2)/8;
    Bs(0,14) = (x + 1)*(-3*pow(x , 2)*y - pow(x , 2) + 3*x*y + x - 15*pow(y,3) - 3*pow(y,2) + 15*y + 3)/8;
    Bs(0,16) = (x + 1)*(y + 1)*(5*pow(y,2) - 2*y - 1)/8;
    Bs(0,17) = (3*y + 1)*(x - 1)*pow(x +1 , 2)/8;
    Bs(0,18) = -(3*y*(x - 1)*(pow(x , 2) + x + 5*pow(y,2) - 5))/8;
    Bs(0,19) = -3*y*pow(x -1 , 2)*(x + 1)/8;
    Bs(0,20) = (x - 1)*(3*pow(x , 2)*y + pow(x , 2) + 3*x*y + x + 15*pow(y,3) + 3*pow(y,2) - 15*y - 3)/8;
    Bs(0,22) = (x - 1)*(y + 1)*(- 5*pow(y,2) + 2*y + 1)/8;
    Bs(0,23) = (3*y + 1)*pow(x -1 , 2)*(x + 1)/8;

    Bs(1,0) = (3*x*(y - 1)*(5*pow(x , 2) + pow(y,2) + y - 5))/8;
    Bs(1,1) = (y - 1)*(15*pow(x,3) - 3*pow(x , 2) + 3*x*pow(y,2) + 3*x*y - 15*x - pow(y,2) - y + 3)/8;
    Bs(1,2) = 3*x*pow((y - 1) , 2)*(y + 1)/8;
    Bs(1,3) = (x - 1)*(y - 1)*(5*pow(x , 2) + 2*x - 1)/8;
    Bs(1,5) = (3*x - 1)*pow((y - 1) , 2)*(y + 1)/8;
    Bs(1,6) = -3*x*(y - 1)*(5*pow(x , 2) + pow(y,2) + y - 5)/8;
    Bs(1,7) = (y - 1)*(15*pow(x,3) + 3*pow(x , 2) + 3*x*pow(y,2) + 3*x*y - 15*x + pow(y,2) + y - 3)/8;
    Bs(1,8)= -3*x*pow((y - 1) , 2)*(y + 1)/8;
    Bs(1,9)= ((x + 1)*(y - 1)*(- 5*pow(x , 2) + 2*x + 1))/8;
    Bs(1,11) = (3*x + 1)*pow((y - 1) , 2)*(y + 1)/8;
    Bs(1,12) = 3*x*(y + 1)*(5*pow(x , 2) + pow(y,2) - y - 5)/8;
    Bs(1,13) = (y + 1)*(-15*pow(x,3) - 3*pow(x , 2) - 3*x*pow(y,2) + 3*x*y + 15*x - pow(y,2) + y + 3)/8;
    Bs(1,14)= -3*x*(y - 1)*pow((y + 1) , 2)/8;
    Bs(1,15) = (x + 1)*(y + 1)*(5*pow(x , 2) - 2*x - 1)/8;
    Bs(1,17) = (3*x + 1)*(y - 1)*pow((y + 1) , 2)/8;
    Bs(1,18) = (3*x*(y + 1)*(-5*pow(x , 2) - pow(y,2) + y + 5))/8;
    Bs(1,19) = (y + 1)*(-15*pow(x,3) + 3*pow(x , 2) - 3*x*pow(y,2) + 3*x*y + 15*x + pow(y,2) - y - 3)/8;
    Bs(1,20) = 3*x*(y - 1)*pow((y + 1) , 2)/8;
    Bs(1,21) = -(x - 1)*(y + 1)*(5*pow(x , 2) + 2*x - 1)/8;
    Bs(1,23) = (3*x - 1)*(y - 1)*pow((y + 1) , 2)/8;

    Bs(2,0) = -(15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bs(2,1) = -(3*x + 1)*(x - 1)*(5*pow(x , 2) + 2*x + 6*pow(y,2) - 9)/32;
    Bs(2,2) = -(3*y + 1)*(y - 1)*(6*pow(x , 2) + 5*pow(y,2) + 2*y - 9)/32;
    Bs(2,3) = -(5*x + 1)*pow(x -1 , 2)*(x + 1)/32;
    Bs(2,4) = -(5*y + 1)*pow((y - 1) , 2)*(y + 1)/32;
    Bs(2,5) = -(3*x + 1)*(3*y + 1)*(x - 1)*(y - 1)/16;
    Bs(2,6) = (15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bs(2,7) = (3*x - 1)*(x + 1)*(- 5*pow(x , 2) + 2*x - 6*pow(y,2) + 9)/32;
    Bs(2,8) = (3*y + 1)*(y - 1)*(6*pow(x , 2) + 5*pow(y,2) + 2*y - 9)/32;
    Bs(2,9) = (5*x - 1)*(x - 1)*pow(x +1 , 2)/32;
    Bs(2,10) = (5*y + 1)*pow((y - 1) , 2)*(y + 1)/32;
    Bs(2,11) = -(3*x - 1)*(3*y + 1)*(x + 1)*(y - 1)/16;
    Bs(2,12) = -(15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bs(2,13) = (3*x - 1)*(x + 1)*(5*pow(x , 2) - 2*x + 6*pow(y,2) - 9)/32;
    Bs(2,14) = (3*y - 1)*(y + 1)*(6*pow(x , 2) + 5*pow(y,2) - 2*y - 9)/32;
    Bs(2,15) = -(5*x - 1)*(x - 1)*pow(x +1 , 2)/32;
    Bs(2,16) = -(5*y - 1)*(y - 1)*pow((y + 1) , 2)/32;
    Bs(2,17) = -(3*x - 1)*(3*y - 1)*(x + 1)*(y + 1)/16;
    Bs(2,18) = (15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bs(2,19) = (3*x + 1)*(x - 1)*(5*pow(x , 2) + 2*x + 6*pow(y,2) - 9)/32;
    Bs(2,20) = (3*y - 1)*(y + 1)*(-6*pow(x , 2) - 5*pow(y,2) + 2*y + 9)/32;
    Bs(2,21) = (5*x + 1)*pow(x -1 , 2)*(x + 1)/32;
    Bs(2,22) = (5*y - 1)*(y - 1)*pow((y + 1) , 2)/32;
    Bs(2,23) = -(3*x + 1)*(3*y - 1)*(x - 1)*(y + 1)/16;


    Eigen::MatrixXd Dmat(3,3);
    double dtdy,dsdy,dtdx,dsdx;
    Eigen::MatrixXd invjac= m_shape.getJacobMatrix().inverse();
    dsdx= invjac(0,0);
    dsdy= invjac(0,1);
    dtdx= invjac(1,0);
    dtdy=invjac(1,1);

    Dmat(0,0)=dtdy*dtdy;
    Dmat(0,1)=dsdy*dsdy;
    Dmat(0,2)=-2*dsdy*dtdy;
    Dmat(1,0)=dtdx*dtdx;
    Dmat(1,1)=dsdx*dsdx;
    Dmat(1,2)=-2*dsdx*dtdx;
    Dmat(2,0)=-dtdx*dtdy;
    Dmat(2,1)=-dsdx*dsdy;
    Dmat(2,2)=dsdx*dtdy+dsdy*dtdx;

    Eigen::MatrixXd J = m_shape.getJacobMatrix();
    double dxds,dxdt,dyds,dydt;
    dxds=J(0,0);dxdt=J(0,1);dyds=J(1,0);dydt=J(1,1);
    Eigen::Matrix<double,6,6> tmat;
    tmat<<  1,0,0,0,0,0,
            0,dxds,dyds,0,0,0,
            0,dxdt,dydt,0,0,0,
            0,0,0,dxds*dxds,dyds*dyds,2*dxds*dyds,
            0,0,0,dxdt*dxdt,dydt*dydt,2*dxdt*dydt,
            0,0,0,dxds*dxdt,dyds*dydt,dxds*dydt+dxdt*dyds;
    Eigen::Matrix<double,6,6>zero;zero.setZero();

    Eigen::Matrix<double,24,24> Tmat;
    Tmat << tmat,zero,zero,zero,
            zero,tmat,zero,zero,
            zero,zero,tmat,zero,
            zero,zero,zero,tmat;
    Oceane::Matrix bsmat(3,24);
    Bs= Dmat*Bs*Tmat;
}


} //namespace Oceane




