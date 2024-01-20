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
#include <vector>

using namespace std;

class PlaneStress: public Oceane::Domain {
    public:
    Oceane::Domain domain;
    vector<double> xi;
    Oceane::nloptsolverplus problem;
    string inputfilename;

    PlaneStress(string filename) {
        inputfilename = filename;
    }

    void read() {
        double temp;
        ifstream input;
        input.open(inputfilename);
        if(input) {
            while(input>>temp)
                xi.push_back(temp);
        }
    } 

    void init_variables() {
        Oceane::AbaqusIO("case00.inp",&domain);
        Oceane::DofManager dofmanager(&domain);

        problem.Init(&dofmanager);
        problem.InitX(xi);

        //creating objective function and adding it to the problem.
        auto objfun = make_shared<Oceane::ObjectiveFunction>(&domain,&dofmanager);
        problem.setObjfun(objfun);

        //creating constraints and adding it to the problem

        for(auto iter: domain.getElements())
        {
            if(iter->isBoundary()) {
                for(auto edge:iter->getEdgeptrs()) {

                    if(edge->isFix()) {
                        auto tx1 =
                                make_shared<Oceane::Boundary::
                                Rectangle::Displacement::Txx>(edge,1);
                        auto tx2 =make_shared<Oceane::Boundary::
                                Rectangle::Displacement::Txx>(edge,2);
                        auto ty1 =make_shared<Oceane::Boundary::
                                Rectangle::Displacement::Tyy>(edge,1);
                        auto ty2 =make_shared<Oceane::Boundary::
                                Rectangle::Displacement::Tyy>(edge,2);
                        auto fx =make_shared<Oceane::Boundary::
                                Rectangle::Displacement::Fxx>(edge);
                        auto fy=make_shared<Oceane::Boundary::
                                Rectangle::Displacement::Fyy>(edge);
                        auto m= make_shared<Oceane::Boundary::
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
                                make_shared<Oceane::Boundary::
                                Rectangle::Traction::Txx>(edge,0);
                        auto tx2 =make_shared<Oceane::Boundary::
                                Rectangle::Traction::Txx>(edge,1);
                        auto ty1 =make_shared<Oceane::Boundary::
                                Rectangle::Traction::Tyy>(edge,0);
                        auto ty2 =make_shared<Oceane::Boundary::
                                Rectangle::Traction::Tyy>(edge,1);
                        auto fx =make_shared<Oceane::Boundary::
                                Rectangle::Traction::Fxx>(edge);
                        auto fy= make_shared<Oceane::Boundary::
                                Rectangle::Traction::Fyy>(edge);
                        auto m= make_shared<Oceane::Boundary::
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
    }

    void writefile(string file, vector<double> vec) {
        ofstream out;
        out.open(file);
        if(out) {
            for(auto& iter:vec)
                out<<setprecision(16)<<iter<<"\n";
        }
    }

    void solve() {
        auto start = chrono::high_resolution_clock::now();
        problem.Solve();
        auto stop =chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout<<"\nRunTime\t"<<duration.count()<<"seconds";
        auto sol =problem.getX();
        auto constraints = problem.getLinearConstraints();

        for(auto& iter: constraints ) {
            double res=0.0;
            iter->Eval(&sol[0],res);
            cout<<"\n"<<res<<"\t"<<iter->getVal();
        }

        cout<<"\n\nSolution";
        for(auto& iter:sol) {
            cout<<"\n"<<iter;
        }
        writefile("res.txt",sol);
    }

};