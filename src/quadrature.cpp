#include "quadrature.h"

namespace Oceane {

Quadrature::Quadrature(int pts)
{
    m_pts =pts;
    m_nodes.resize(m_pts);
    m_weight.resize(m_pts);
    switch(m_pts)
    {
    case 2:
    {
        double x2[1] = {0.5773502691896257645091488};
        double w2[1] = {1.0000000000000000000000000};
        set_gauss_val(x2,w2);
        break;
    }
    case 4:
    {
        double x4[2] = {0.3399810435848562648026658,0.8611363115940525752239465};
        double w4[2] = {0.6521451548625461426269361,0.3478548451374538573730639};
        set_gauss_val(x4,w4);
        break;
    }
    case 6:
    {
        double x6[3] = {0.2386191860831969086305017,0.6612093864662645136613996,0.9324695142031520278123016};
        double w6[3] = {0.4679139345726910473898703,0.3607615730481386075698335,0.1713244923791703450402961};
        set_gauss_val(x6,w6);
        break;
    }
    case 8:
    {
        double x8[4] = {0.1834346424956498049394761,0.5255324099163289858177390,0.7966664774136267395915539,0.9602898564975362316835609};
        double w8[4] = {0.3626837833783619829651504,0.3137066458778872873379622,0.2223810344533744705443560,0.1012285362903762591525314};
        set_gauss_val(x8,w8);
        break;
    }
    case 10:
    {
        double x10[5] = {0.1488743389816312108848260,0.4333953941292471907992659,0.6794095682990244062343274,0.8650633666889845107320967,0.9739065285171717200779640};
        double w10[5] = {0.2955242247147528701738930,0.2692667193099963550912269,0.2190863625159820439955349,0.1494513491505805931457763,0.0666713443086881375935688};
        set_gauss_val(x10,w10);
        break;
    }
    case 12:
    {
        double x12[6] = {0.1252334085114689154724414,0.3678314989981801937526915,0.5873179542866174472967024,0.7699026741943046870368938,0.9041172563704748566784659,0.9815606342467192506905491};
        double w12[6] = {0.2491470458134027850005624,0.2334925365383548087608499,0.2031674267230659217490645,0.1600783285433462263346525,0.1069393259953184309602547,0.0471753363865118271946160};
        set_gauss_val(x12,w12);
        break;
    }
    case 14:
    {
        double x14[7] = {0.1080549487073436620662447,0.3191123689278897604356718,0.5152486363581540919652907,0.6872929048116854701480198,0.8272013150697649931897947,0.9284348836635735173363911,0.9862838086968123388415973};
        double w14[7] = {0.2152638534631577901958764,0.2051984637212956039659241,0.1855383974779378137417166,0.1572031671581935345696019,0.1215185706879031846894148,0.0801580871597602098056333,0.0351194603317518630318329};
        set_gauss_val(x14,w14);
        break;
    }

    default:
        std::cerr<<"This quadrature rule isn't implemented yet";
    }
}



void Quadrature::set_gauss_val(const double* x,const double* w)
{
    for(int i=0; i<m_pts/2; i++)
    {
        m_nodes[2*i]=-1.0*x[i];
        m_nodes[2*i+1]=x[i];
        m_weight[2*i]=1.0*w[i];
        m_weight[2*i+1]=w[i];

       /* auto val=-1.0*x[i];
        m_nodes.push_back(val);
        m_nodes.push_back(x[i]);
        m_weight.push_back(-1.0*w[i]);
        m_weight.push_back(w[i]);*/
    }
}

double Quadrature::AreaIntegral( double (*function)(double, double))
{
    double val=0;
    for(int i=0;i<m_nodes.size();i++)
        for(int j=0;j<m_nodes.size();j++)
        {
            val+=m_weight[i]*m_weight[j]*function(m_nodes[i],m_nodes[j]);
        }
    return val;
}
/*
Eigen::MatrixXd Quadrature::AreaIntegral(int, Eigen::MatrixXd (*function)(double, double))
{
    Eigen::MatrixXd val;
    val.setZero();
    for(int i=0;i<m_nodes.size();i++)
        for(int j=0;j<m_nodes.size();j++)
        {
            val+=m_weight[i]*m_weight[j]*function(m_nodes[i],m_nodes[j]);
        }
    return val;
}
*/


} // namespace Oceane
