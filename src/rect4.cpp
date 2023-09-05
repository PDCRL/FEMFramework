#include "rect4.h"

namespace Oceane {

Rect4::Rect4(Index id,std::vector<Nodeptr> nodes)
    :Element (id,nodes)
{
    Oceane::Vector xvec(4);
    Oceane::Vector yvec(4);
    int xi=-1,yi=-1;
    for(auto node:m_nodes)
    {
        auto coord =node->getCoord();
        xvec[++xi]=coord[0];
        yvec[++yi]=coord[1];
    }
    m_xcoord =xvec;
    m_ycoord =yvec;
}


void Rect4::getBsmat_at(double x, double y, Matrix& Bs, double &detJ)
{
    Oceane::Matrix JacobMatrix(2,2);
    Oceane::Vector dNx(4);
    Oceane::Vector dNy(4);


    dNx[0]=(x-1)/4;
    dNx[1]=-(x-1)/4;
    dNx[2]=(x+1)/4;
    dNx[3]=-(x+1)/4;

    dNy[0]= (y-1)/4;
    dNy[1]=-(y+1)/4;
    dNy[2]=(y+1)/4;
    dNy[3]=-(y-1)/4;

    JacobMatrix(0,0) = dNx.transpose() * m_xcoord;
    JacobMatrix(0,1) = dNy.transpose() * m_xcoord;
    JacobMatrix(1,0) = dNx.transpose() * m_ycoord;
    JacobMatrix(1,1) = dNy.transpose() * m_ycoord;

    detJ = JacobMatrix.determinant();

    if(Bs.rows()!=3 or Bs.cols()!=24)
        Bs.resize(3,24);

    Oceane::Matrix Bsmat(3,24);
    Bsmat.setZero();

    Bsmat.setZero();
    Bsmat(0,0)=(3*y*(x - 1)*(pow(x,2) + x + 5*pow(y,2) - 5))/8;
    Bsmat(0,1)= 3*y*pow(x -1 , 2)*(x + 1)/8;
    Bsmat(0,2) = (x - 1)*(3*pow(x , 2)*y - pow(x , 2) + 3*x*y - x + 15*pow(y,3) - 3*pow(y,2) - 15*y + 3)/8;
    Bsmat(0,4) = (x - 1)*(y - 1)*(5*pow(y,2) + 2*y - 1)/8;
    Bsmat(0,5) = (3*y - 1)*pow(x -1 , 2)*(x + 1)/8;
    Bsmat(0,6) = 3*y*(x + 1)*(- pow(x , 2) + x - 5*pow(y,2) + 5)/8;
    Bsmat(0,7) = 3*y*(x - 1)*pow(x +1 , 2)/8;
    Bsmat(0,8) = -(x + 1)*(3*pow(x , 2)*y - pow(x , 2) - 3*x*y + x + 15*pow(y,3) - 3*pow(y,2) - 15*y + 3)/8;
    Bsmat(0,10) = -(x + 1)*(y - 1)*(5*pow(y,2) + 2*y - 1)/8;
    Bsmat(0,11) = (3*y - 1)*(x - 1)*pow(x +1 , 2)/8;
    Bsmat(0,12) = 3*y*(x + 1)*(pow(x , 2) - x + 5*pow(y,2) - 5)/8;
    Bsmat(0,13) = -3*y*(x - 1)*pow(x +1 , 2)/8;
    Bsmat(0,14) = (x + 1)*(-3*pow(x , 2)*y - pow(x , 2) + 3*x*y + x - 15*pow(y,3) - 3*pow(y,2) + 15*y + 3)/8;
    Bsmat(0,16) = (x + 1)*(y + 1)*(5*pow(y,2) - 2*y - 1)/8;
    Bsmat(0,17) = (3*y + 1)*(x - 1)*pow(x +1 , 2)/8;
    Bsmat(0,18) = -(3*y*(x - 1)*(pow(x , 2) + x + 5*pow(y,2) - 5))/8;
    Bsmat(0,19) = -3*y*pow(x -1 , 2)*(x + 1)/8;
    Bsmat(0,20) = (x - 1)*(3*pow(x , 2)*y + pow(x , 2) + 3*x*y + x + 15*pow(y,3) + 3*pow(y,2) - 15*y - 3)/8;
    Bsmat(0,22) = (x - 1)*(y + 1)*(- 5*pow(y,2) + 2*y + 1)/8;
    Bsmat(0,23) = (3*y + 1)*pow(x -1 , 2)*(x + 1)/8;

    Bsmat(1,0) = (3*x*(y - 1)*(5*pow(x , 2) + pow(y,2) + y - 5))/8;
    Bsmat(1,1) = (y - 1)*(15*pow(x,3) - 3*pow(x , 2) + 3*x*pow(y,2) + 3*x*y - 15*x - pow(y,2) - y + 3)/8;
    Bsmat(1,2) = 3*x*pow((y - 1) , 2)*(y + 1)/8;
    Bsmat(1,3) = (x - 1)*(y - 1)*(5*pow(x , 2) + 2*x - 1)/8;
    Bsmat(1,5) = (3*x - 1)*pow((y - 1) , 2)*(y + 1)/8;
    Bsmat(1,6) = -3*x*(y - 1)*(5*pow(x , 2) + pow(y,2) + y - 5)/8;
    Bsmat(1,7) = (y - 1)*(15*pow(x,3) + 3*pow(x , 2) + 3*x*pow(y,2) + 3*x*y - 15*x + pow(y,2) + y - 3)/8;
    Bsmat(1,8)= -3*x*pow((y - 1) , 2)*(y + 1)/8;
    Bsmat(1,9)= ((x + 1)*(y - 1)*(- 5*pow(x , 2) + 2*x + 1))/8;
    Bsmat(1,11) = (3*x + 1)*pow((y - 1) , 2)*(y + 1)/8;
    Bsmat(1,12) = 3*x*(y + 1)*(5*pow(x , 2) + pow(y,2) - y - 5)/8;
    Bsmat(1,13) = (y + 1)*(-15*pow(x,3) - 3*pow(x , 2) - 3*x*pow(y,2) + 3*x*y + 15*x - pow(y,2) + y + 3)/8;
    Bsmat(1,14)= -3*x*(y - 1)*pow((y + 1) , 2)/8;
    Bsmat(1,15) = (x + 1)*(y + 1)*(5*pow(x , 2) - 2*x - 1)/8;
    Bsmat(1,17) = (3*x + 1)*(y - 1)*pow((y + 1) , 2)/8;
    Bsmat(1,18) = (3*x*(y + 1)*(-5*pow(x , 2) - pow(y,2) + y + 5))/8;
    Bsmat(1,19) = (y + 1)*(-15*pow(x,3) + 3*pow(x , 2) - 3*x*pow(y,2) + 3*x*y + 15*x + pow(y,2) - y - 3)/8;
    Bsmat(1,20) = 3*x*(y - 1)*pow((y + 1) , 2)/8;
    Bsmat(1,21) = -(x - 1)*(y + 1)*(5*pow(x , 2) + 2*x - 1)/8;
    Bsmat(1,23) = (3*x - 1)*(y - 1)*pow((y + 1) , 2)/8;

    Bsmat(2,0) = -(15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bsmat(2,1) = -(3*x + 1)*(x - 1)*(5*pow(x , 2) + 2*x + 6*pow(y,2) - 9)/32;
    Bsmat(2,2) = -(3*y + 1)*(y - 1)*(6*pow(x , 2) + 5*pow(y,2) + 2*y - 9)/32;
    Bsmat(2,3) = -(5*x + 1)*pow(x -1 , 2)*(x + 1)/32;
    Bsmat(2,4) = -(5*y + 1)*pow((y - 1) , 2)*(y + 1)/32;
    Bsmat(2,5) = -(3*x + 1)*(3*y + 1)*(x - 1)*(y - 1)/16;
    Bsmat(2,6) = (15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bsmat(2,7) = (3*x - 1)*(x + 1)*(- 5*pow(x , 2) + 2*x - 6*pow(y,2) + 9)/32;
    Bsmat(2,8) = (3*y + 1)*(y - 1)*(6*pow(x , 2) + 5*pow(y,2) + 2*y - 9)/32;
    Bsmat(2,9) = (5*x - 1)*(x - 1)*pow(x +1 , 2)/32;
    Bsmat(2,10) = (5*y + 1)*pow((y - 1) , 2)*(y + 1)/32;
    Bsmat(2,11) = -(3*x - 1)*(3*y + 1)*(x + 1)*(y - 1)/16;
    Bsmat(2,12) = -(15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bsmat(2,13) = (3*x - 1)*(x + 1)*(5*pow(x , 2) - 2*x + 6*pow(y,2) - 9)/32;
    Bsmat(2,14) = (3*y - 1)*(y + 1)*(6*pow(x , 2) + 5*pow(y,2) - 2*y - 9)/32;
    Bsmat(2,15) = -(5*x - 1)*(x - 1)*pow(x +1 , 2)/32;
    Bsmat(2,16) = -(5*y - 1)*(y - 1)*pow((y + 1) , 2)/32;
    Bsmat(2,17) = -(3*x - 1)*(3*y - 1)*(x + 1)*(y + 1)/16;
    Bsmat(2,18) = (15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bsmat(2,19) = (3*x + 1)*(x - 1)*(5*pow(x , 2) + 2*x + 6*pow(y,2) - 9)/32;
    Bsmat(2,20) = (3*y - 1)*(y + 1)*(-6*pow(x , 2) - 5*pow(y,2) + 2*y + 9)/32;
    Bsmat(2,21) = (5*x + 1)*pow(x -1 , 2)*(x + 1)/32;
    Bsmat(2,22) = (5*y - 1)*(y - 1)*pow((y + 1) , 2)/32;
    Bsmat(2,23) = -(3*x + 1)*(3*y - 1)*(x - 1)*(y + 1)/16;





    Eigen::MatrixXd Dmat(3,3);
    double dtdy,dsdy,dtdx,dsdx;
    Eigen::MatrixXd invjac= JacobMatrix.inverse();
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




    Eigen::MatrixXd J = JacobMatrix;
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
    Bs= Dmat*Bsmat*Tmat;

}

void Rect4::getBdmat_at(double x, double y, Matrix &Bdmat, double &detJ)
{

    Oceane::Matrix JacobMatrix(2,2);
    Oceane::Vector dNx(4);
    Oceane::Vector dNy(4);


    dNx[0]=(x-1)/4;
    dNx[1]=-(x-1)/4;
    dNx[2]=(x+1)/4;
    dNx[3]=-(x+1)/4;

    dNy[0]= (y-1)/4;
    dNy[1]=-(y+1)/4;
    dNy[2]=(y+1)/4;
    dNy[3]=-(y-1)/4;

    JacobMatrix(0,0) = dNx.transpose() * m_xcoord;
    JacobMatrix(0,1) = dNy.transpose() * m_xcoord;
    JacobMatrix(1,0) = dNx.transpose() * m_ycoord;
    JacobMatrix(1,1) = dNy.transpose() * m_ycoord;

    detJ = JacobMatrix.determinant();

    if(Bdmat.rows()!=3 or (Bdmat.cols()!=8)) Bdmat.resize(3,8);

    Eigen::MatrixXd inv= JacobMatrix.inverse();
    Eigen::VectorXd Nis,Nit;
    Nis=dNx;
    Nit=dNy;
    Bdmat.resize(3,8);
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

void Rect4::getBmats_at(double x, double y, Matrix& Bdmat,Matrix &Bs, double &detJ)
{
    Oceane::Matrix JacobMatrix(2,2);
    Oceane::Vector dNx(4);
    Oceane::Vector dNy(4);


    dNx[0]=(x-1)/4;
    dNx[1]=-(x-1)/4;
    dNx[2]=(x+1)/4;
    dNx[3]=-(x+1)/4;

    dNy[0]= (y-1)/4;
    dNy[1]=-(y+1)/4;
    dNy[2]=(y+1)/4;
    dNy[3]=-(y-1)/4;

    JacobMatrix(0,0) = dNx.transpose() * m_xcoord;
    JacobMatrix(0,1) = dNy.transpose() * m_xcoord;
    JacobMatrix(1,0) = dNx.transpose() * m_ycoord;
    JacobMatrix(1,1) = dNy.transpose() * m_ycoord;

    detJ = JacobMatrix.determinant();

    if(Bs.rows()!=3 or Bs.cols()!=24)
        Bs.resize(3,24);

    Oceane::Matrix Bsmat(3,24);
    Bsmat.setZero();

    Bsmat.setZero();
    Bsmat(0,0)=(3*y*(x - 1)*(pow(x,2) + x + 5*pow(y,2) - 5))/8;
    Bsmat(0,1)= 3*y*pow(x -1 , 2)*(x + 1)/8;
    Bsmat(0,2) = (x - 1)*(3*pow(x , 2)*y - pow(x , 2) + 3*x*y - x + 15*pow(y,3) - 3*pow(y,2) - 15*y + 3)/8;
    Bsmat(0,4) = (x - 1)*(y - 1)*(5*pow(y,2) + 2*y - 1)/8;
    Bsmat(0,5) = (3*y - 1)*pow(x -1 , 2)*(x + 1)/8;
    Bsmat(0,6) = 3*y*(x + 1)*(- pow(x , 2) + x - 5*pow(y,2) + 5)/8;
    Bsmat(0,7) = 3*y*(x - 1)*pow(x +1 , 2)/8;
    Bsmat(0,8) = -(x + 1)*(3*pow(x , 2)*y - pow(x , 2) - 3*x*y + x + 15*pow(y,3) - 3*pow(y,2) - 15*y + 3)/8;
    Bsmat(0,10) = -(x + 1)*(y - 1)*(5*pow(y,2) + 2*y - 1)/8;
    Bsmat(0,11) = (3*y - 1)*(x - 1)*pow(x +1 , 2)/8;
    Bsmat(0,12) = 3*y*(x + 1)*(pow(x , 2) - x + 5*pow(y,2) - 5)/8;
    Bsmat(0,13) = -3*y*(x - 1)*pow(x +1 , 2)/8;
    Bsmat(0,14) = (x + 1)*(-3*pow(x , 2)*y - pow(x , 2) + 3*x*y + x - 15*pow(y,3) - 3*pow(y,2) + 15*y + 3)/8;
    Bsmat(0,16) = (x + 1)*(y + 1)*(5*pow(y,2) - 2*y - 1)/8;
    Bsmat(0,17) = (3*y + 1)*(x - 1)*pow(x +1 , 2)/8;
    Bsmat(0,18) = -(3*y*(x - 1)*(pow(x , 2) + x + 5*pow(y,2) - 5))/8;
    Bsmat(0,19) = -3*y*pow(x -1 , 2)*(x + 1)/8;
    Bsmat(0,20) = (x - 1)*(3*pow(x , 2)*y + pow(x , 2) + 3*x*y + x + 15*pow(y,3) + 3*pow(y,2) - 15*y - 3)/8;
    Bsmat(0,22) = (x - 1)*(y + 1)*(- 5*pow(y,2) + 2*y + 1)/8;
    Bsmat(0,23) = (3*y + 1)*pow(x -1 , 2)*(x + 1)/8;

    Bsmat(1,0) = (3*x*(y - 1)*(5*pow(x , 2) + pow(y,2) + y - 5))/8;
    Bsmat(1,1) = (y - 1)*(15*pow(x,3) - 3*pow(x , 2) + 3*x*pow(y,2) + 3*x*y - 15*x - pow(y,2) - y + 3)/8;
    Bsmat(1,2) = 3*x*pow((y - 1) , 2)*(y + 1)/8;
    Bsmat(1,3) = (x - 1)*(y - 1)*(5*pow(x , 2) + 2*x - 1)/8;
    Bsmat(1,5) = (3*x - 1)*pow((y - 1) , 2)*(y + 1)/8;
    Bsmat(1,6) = -3*x*(y - 1)*(5*pow(x , 2) + pow(y,2) + y - 5)/8;
    Bsmat(1,7) = (y - 1)*(15*pow(x,3) + 3*pow(x , 2) + 3*x*pow(y,2) + 3*x*y - 15*x + pow(y,2) + y - 3)/8;
    Bsmat(1,8)= -3*x*pow((y - 1) , 2)*(y + 1)/8;
    Bsmat(1,9)= ((x + 1)*(y - 1)*(- 5*pow(x , 2) + 2*x + 1))/8;
    Bsmat(1,11) = (3*x + 1)*pow((y - 1) , 2)*(y + 1)/8;
    Bsmat(1,12) = 3*x*(y + 1)*(5*pow(x , 2) + pow(y,2) - y - 5)/8;
    Bsmat(1,13) = (y + 1)*(-15*pow(x,3) - 3*pow(x , 2) - 3*x*pow(y,2) + 3*x*y + 15*x - pow(y,2) + y + 3)/8;
    Bsmat(1,14)= -3*x*(y - 1)*pow((y + 1) , 2)/8;
    Bsmat(1,15) = (x + 1)*(y + 1)*(5*pow(x , 2) - 2*x - 1)/8;
    Bsmat(1,17) = (3*x + 1)*(y - 1)*pow((y + 1) , 2)/8;
    Bsmat(1,18) = (3*x*(y + 1)*(-5*pow(x , 2) - pow(y,2) + y + 5))/8;
    Bsmat(1,19) = (y + 1)*(-15*pow(x,3) + 3*pow(x , 2) - 3*x*pow(y,2) + 3*x*y + 15*x + pow(y,2) - y - 3)/8;
    Bsmat(1,20) = 3*x*(y - 1)*pow((y + 1) , 2)/8;
    Bsmat(1,21) = -(x - 1)*(y + 1)*(5*pow(x , 2) + 2*x - 1)/8;
    Bsmat(1,23) = (3*x - 1)*(y - 1)*pow((y + 1) , 2)/8;

    Bsmat(2,0) = -(15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bsmat(2,1) = -(3*x + 1)*(x - 1)*(5*pow(x , 2) + 2*x + 6*pow(y,2) - 9)/32;
    Bsmat(2,2) = -(3*y + 1)*(y - 1)*(6*pow(x , 2) + 5*pow(y,2) + 2*y - 9)/32;
    Bsmat(2,3) = -(5*x + 1)*pow(x -1 , 2)*(x + 1)/32;
    Bsmat(2,4) = -(5*y + 1)*pow((y - 1) , 2)*(y + 1)/32;
    Bsmat(2,5) = -(3*x + 1)*(3*y + 1)*(x - 1)*(y - 1)/16;
    Bsmat(2,6) = (15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bsmat(2,7) = (3*x - 1)*(x + 1)*(- 5*pow(x , 2) + 2*x - 6*pow(y,2) + 9)/32;
    Bsmat(2,8) = (3*y + 1)*(y - 1)*(6*pow(x , 2) + 5*pow(y,2) + 2*y - 9)/32;
    Bsmat(2,9) = (5*x - 1)*(x - 1)*pow(x +1 , 2)/32;
    Bsmat(2,10) = (5*y + 1)*pow((y - 1) , 2)*(y + 1)/32;
    Bsmat(2,11) = -(3*x - 1)*(3*y + 1)*(x + 1)*(y - 1)/16;
    Bsmat(2,12) = -(15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bsmat(2,13) = (3*x - 1)*(x + 1)*(5*pow(x , 2) - 2*x + 6*pow(y,2) - 9)/32;
    Bsmat(2,14) = (3*y - 1)*(y + 1)*(6*pow(x , 2) + 5*pow(y,2) - 2*y - 9)/32;
    Bsmat(2,15) = -(5*x - 1)*(x - 1)*pow(x +1 , 2)/32;
    Bsmat(2,16) = -(5*y - 1)*(y - 1)*pow((y + 1) , 2)/32;
    Bsmat(2,17) = -(3*x - 1)*(3*y - 1)*(x + 1)*(y + 1)/16;
    Bsmat(2,18) = (15*pow(x , 4) + 18*pow(x , 2)*pow(y,2) - 36*pow(x , 2) + 15*pow(y , 4) - 36*pow(y,2) + 24)/32;
    Bsmat(2,19) = (3*x + 1)*(x - 1)*(5*pow(x , 2) + 2*x + 6*pow(y,2) - 9)/32;
    Bsmat(2,20) = (3*y - 1)*(y + 1)*(-6*pow(x , 2) - 5*pow(y,2) + 2*y + 9)/32;
    Bsmat(2,21) = (5*x + 1)*pow(x -1 , 2)*(x + 1)/32;
    Bsmat(2,22) = (5*y - 1)*(y - 1)*pow((y + 1) , 2)/32;
    Bsmat(2,23) = -(3*x + 1)*(3*y - 1)*(x - 1)*(y + 1)/16;





    Eigen::MatrixXd Dmat(3,3);
    double dtdy,dsdy,dtdx,dsdx;
    Eigen::MatrixXd invjac= JacobMatrix.inverse();
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




    Eigen::MatrixXd J = JacobMatrix;
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
    Bs= Dmat*Bsmat*Tmat;

    if(Bdmat.rows()!=3 or (Bdmat.cols()!=8)) Bdmat.resize(3,8);

    Eigen::MatrixXd inv= JacobMatrix.inverse();
    Eigen::VectorXd Nis,Nit;
    Nis=dNx;
    Nit=dNy;
    Bdmat.resize(3,8);
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

} // namespace Oceane
