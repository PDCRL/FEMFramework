#include <iostream>
#include "domain.h"
#include "abaqusio.h"
#include "dofmanager.h"
#include "objectivefunction.h"
#include "rectboundary.h"
#include "nloptsolver.h"
#include "nloptsolverplus.h"
#include <ctime>
#include <iostream>
#include <chrono>
#include <thread>
#include <assert.h>
#include <iomanip>

std::vector<double> readfile(std::string file)
{
    std::vector<double> xi; double temp;
    std::ifstream input;
    input.open(file);
    if(input) {
        while(input>>temp)
            xi.push_back(temp);
    }
    return xi;
}

void writefile(std::string file, std::vector<double> vec)
{
    std::ofstream out;
    out.open(file);
    if(out) {
        for(auto& iter:vec)
            out<<std::setprecision(16)<<iter<<"\n";
    }
}
int main()
{
    auto xi = readfile("test1.txt");

    Oceane::Domain domain;
    Oceane::AbaqusIO("case00.inp",&domain);
    Oceane::DofManager dofmanager(&domain);

//    Oceane::WorhpSolver problem;
    Oceane::nloptsolverplus problem;
//    Oceane::nloptsolver problem;
//    Oceane::IpoptSolver problem;
    problem.Init(&dofmanager);
    problem.Ini tX(xi);

    //creating objective function and adding it to the problem.
//    auto objfun = std::make_shared<Oceane::SquareNorm>(&domain,&dofmanager);
    auto objfun = std::make_shared<Oceane::ObjectiveFunction>(&domain,&dofmanager);
    problem.setObjfun(objfun);


    //creating constraints and adding it to the problem

    for(auto iter: domain.getElements())
    {
        if(iter->isBoundary()) {
            for(auto edge:iter->getEdgeptrs()) {

                if(edge->isFix()) {
                    auto tx1 =
                            std::make_shared<Oceane::Boundary::
                            Rectangle::Displacement::Txx>(edge,1);
                    auto tx2 =std::make_shared<Oceane::Boundary::
                            Rectangle::Displacement::Txx>(edge,2);
                    auto ty1 =std::make_shared<Oceane::Boundary::
                            Rectangle::Displacement::Tyy>(edge,1);
                    auto ty2 =std::make_shared<Oceane::Boundary::
                            Rectangle::Displacement::Tyy>(edge,2);
                    auto fx =std::make_shared<Oceane::Boundary::
                            Rectangle::Displacement::Fxx>(edge);
                    auto fy= std::make_shared<Oceane::Boundary::
                            Rectangle::Displacement::Fyy>(edge);
                    auto m= std::make_shared<Oceane::Boundary::
                            Rectangle::Displacement::Moment>(edge);
                    problem.addConstraint(tx1);
                    problem.addConstraint(tx2);
                    problem.addConstraint(ty1);
                    problem.addConstraint(ty2);
                    problem.addConstraint(fx);
                    problem.addConstraint(m);
                }
                else {
                    auto tx1 =
                            std::make_shared<Oceane::Boundary::
                            Rectangle::Traction::Txx>(edge,0);
                    auto tx2 =std::make_shared<Oceane::Boundary::
                            Rectangle::Traction::Txx>(edge,1);
                    auto ty1 =std::make_shared<Oceane::Boundary::
                            Rectangle::Traction::Tyy>(edge,0);
                    auto ty2 =std::make_shared<Oceane::Boundary::
                            Rectangle::Traction::Tyy>(edge,1);
                    auto fx =std::make_shared<Oceane::Boundary::
                            Rectangle::Traction::Fxx>(edge);
                    auto fy= std::make_shared<Oceane::Boundary::
                            Rectangle::Traction::Fyy>(edge);
                    auto m= std::make_shared<Oceane::Boundary::
                            Rectangle::Traction::Moment>(edge);
                    problem.addConstraint(tx1);
                    problem.addConstraint(tx2);
                    problem.addConstraint(ty1);
                    problem.addConstraint(ty2);
                    problem.addConstraint(fx);
                    problem.addConstraint(fy);
                    problem.addConstraint(m);
                }
            }
        }
    }


    auto start = std::chrono::high_resolution_clock::now();
    problem.Solve();
    auto stop =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout<<"\nRunTime\t"<<duration.count()<<"seconds";
    auto sol =problem.getX();
    auto constraints = problem.getLinearConstraints();
   for(auto& iter: constraints )
   {
       double res=0.0;
       iter->Eval(&sol[0],res);
       std::cout<<"\n"<<res<<"\t"<<iter->getVal();
   }
   std::cout<<"\n\nSolution";
   for(auto& iter:sol)
   {
       std::cout<<"\n"<<iter;

   }

       writefile("res.txt",sol);

}
