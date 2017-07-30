
 #include <iostream>
#include <stdexcept>
#include "capd/capdlib.h"
#include "capd/krak/krak.h"
#include "capd/mpcapdlib.h"
using namespace std;
using namespace capd;
using capd::autodiff::Node;




void getVbad(Node t, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node
params[], int /*noParams*/)
{
  Node ce = params[0];
  Node cs = params[1];

  Node xSquare = in[0]*in[0];
  Node ySquare = in[1]*in[1];
  Node zSquare = in[2]*in[2];

  Node s = sqrt((xSquare + ySquare + zSquare)*(xSquare + ySquare + zSquare)*(xSquare
+ ySquare + zSquare));
  Node trzy = Node(3.);
  Node jeden = Node(1.);

  out[3] = in[0]*(trzy*cs*(in[5]*in[1]-in[4]*in[2])-ce)/s;
  out[4] = (in[5]*(cs-trzy*cs*xSquare)    - ce*in[1] + trzy*cs*in[3]*in[0]*in[2])/s;
  out[5] = (cs*in[4]*(trzy*xSquare-jeden) - ce*in[2] - trzy*cs*in[3]*in[0]*in[1])/s;

  out[0] = in[3]+Node(0.);
  out[1] = in[4]+Node(0.);
  out[2] = in[5]+Node(0.);
}

void getV(Node t, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node
params[], int /*noParams*/)
{
    Node x = in[0];
    Node y = in[1];
    Node z = in[2];

    Node vx = in[3];
    Node vy = in[4];
    Node vz = in[5];

    Node ce = params[0];
    Node cs = params[1];

    Node minuss = Node(-1);

    Node xSquare = x^2;
    Node ySquare = y^2;
    Node zSquare = z^2;

    Node squaresSum = xSquare + ySquare + zSquare;
    Node con = Node(0.000000979814);

    Node partIIx = cs * (((minuss*3*z)*(minuss*vz*y + vy*z))*((squaresSum)^-2.5) + (vy*((squaresSum)^-1.5)));
    Node partIIy = cs * (((minuss*3*z)*(vz*x - vx*z))*((squaresSum)^-2.5) - (vx*((squaresSum)^-1.5)));
    Node partIIz = Node(0);//(ce*z)/((squaresSum)^1.5);

    out[0] = vx;
    out[1] = vy;
    out[2] = vz;

    out[3] = minuss * ((ce * x) * (squaresSum^-1.5 )  ) - partIIx;
    out[4] = minuss * ((ce * y) * (squaresSum^-1.5 )  ) - partIIy;
    out[5] = (con *(minuss*vy*x + vx*y)*(z))*((squaresSum)^-2.5) - partIIz;
}

void myGetV3(Node t, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node
params[], int /*noParams*/){

    Node squaresSum3 = (sqr(in[0]) + sqr(in[1]) + sqr(in[2])) ^ 1.5;
    Node squaresSum5 = (sqr(in[0]) + sqr(in[1]) + sqr(in[2])) ^ 2.5;

    out[0] = in[3];
    out[1] = in[4];
    out[2] = in[5];

    out[3] = (-in[0] * params[0] ) / squaresSum3 + (Node(6.) * params[1] * in[2] * (Node(-1.) * in[5] * in[1] + in[4] * in[2])) / (squaresSum5) - Node(2.) * params[1] * in[4]  / squaresSum3;
    out[4] = (-in[1] * params[0] ) / squaresSum3 + (Node(6.) * params[1] * in[2] * ( in[5] * in[0] - in[3] * in[2])) / (squaresSum5) + Node(2.) * params[1] * in[3]  / squaresSum3;
    out[5] = Node(6.) * params[1] * (Node(-1.) * in[4] * in[0] + in[3] * in[1]) * in[2]/ squaresSum5 - params[0] * in[2] / squaresSum3;

}

void myGetV(Node t, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node
params[], int /*noParams*/){
    Node x = in[0];
    Node y = in[1];
    Node z = in[2];

    Node vx = in[3];
    Node vy = in[4];
    Node vz = in[5];

    Node ce = params[0];
    Node cs = params[1];


    Node xSquare = sqr(x);
    Node ySquare = sqr(y);
    Node zSquare = sqr(z);

    //Node squaresSum = xSquare + ySquare + zSquare;
    Node squaresSum3 = (sqr(in[0]) + sqr(in[1]) + sqr(in[2])) ^ 1.5;
    Node squaresSum5 = (sqr(in[0]) + sqr(in[1]) + sqr(in[2])) ^ 2.5;

    //Node s3 = (squaresSum^1.5);
    //Node s5 = ()squaresSum^5);

    out[0] = vx;
    out[1] = vy;
    out[2] = Node(0.);

    //out[4] = (Node(-1.) * ce * x / s3) - Node(2.) * cs * (vy / s3);
    //out[5] = (Node(-1.) * ce * y / s3) + Node(2.) * cs * (vx / s3);
    out[3] = (-in[0] * params[0] + Node(-2.) * params[1] * in[4] ) / squaresSum3;
    out[4] = (-in[1] * params[0] - Node(-2.) * params[1] * in[3] ) / squaresSum3;
    //out[5] = (con *(minuss*vy*x + vx*y)*(z))*((squaresSum)^-2.5) - partIIz;
    out[5] = Node(0.);


}

void getV2(Node t, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node squaresSum = (sqr(in[0]) + sqr(in[1]))^1.5;
  out[0] = in[2];
  out[1] = in[3];
  out[2] = (-in[0] * params[0] + Node(-2.) * params[1] * in[3] ) / squaresSum;
  out[3] = (-in[1] * params[0] - Node(-2.) * params[1] * in[2] ) / squaresSum;
}

template<class T, class S, class F>
void show(T t, const S& s, F& f){
  auto u = s(t);
  cout << "t: " << t << " " << u << " " << f(u)<< endl;
}

namespace LD{
  typedef long double Real;
  typedef LDMap Map;
  typedef LDOdeSolver OdeSolver;
  typedef LDTimeMap TimeMap;
  typedef LDVector Vector;
}

namespace Mp{
  typedef MpFloat Real;
  typedef MpMap Map;
  typedef MpOdeSolver OdeSolver;
  typedef MpTimeMap TimeMap;
  typedef MpVector Vector;
}

int main(int, char*[]){

  using namespace LD;
  MpFloat::setDefaultPrecision(1280);
  cout.precision(10);
  int paramsNumber = 2;
  int dim = 6;

  // funkcja getV2 jest identyczna jak getV tylko jest uproszczona - os z nie jest uwzgledniana
  Map f(dim==4? getV2 : myGetV3, dim, dim, paramsNumber);

  //Map f(getV,dim,dim,paramsNumber);

  Real CS = 1.633022997504943e-7;
  Real CE = 0.25326384620828446;

  f.setParameter(0,CE);
  f.setParameter(1,CS);

  int order = 20;
  OdeSolver s(f,order, OdeSolver::StepControlType(2,1e-12));
  s.setAbsoluteTolerance(1e-20);
  s.setRelativeTolerance(1e-20);

  Real initTime = 0.0;
  Real finalTime = 30.5;
  Real startTime = 0.0;

  TimeMap::SolutionCurve solution(initTime);
  TimeMap tm3(s);
  Real d[] = {Real(1.0),Real(0.0),Real(0.0),Real(0.0),Real(0.0),Real(0.0)};


  Vector u3(dim,d);
  cout << f(u3) << endl;

  tm3(finalTime, u3,solution);
  cout << "domain = [" << solution.getLeftDomain() << "," <<
  solution.getRightDomain() << "]\n";
  cout << "------" << endl;

  int N = 500;

//  for(int i=0;i<=N;++i){
  //show(i*solution.getRightDomain()/N,solution,f);
  //}


    for(double d = 0.0 ; d < 30.5; d = d + 0.1){
        cout << "t: " << d << " : "<< solution(d) << endl;
    }

  cout << "------" << endl;
}
