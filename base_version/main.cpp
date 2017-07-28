#include <iostream>
#include <stdexcept>
#include "capd/capdlib.h"
#include "capd/krak/krak.h"
#include "capd/mpcapdlib.h"
using namespace std;
using namespace capd;
using capd::autodiff::Node;




void getV(Node t, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node ce = params[0];
  Node cs = params[1];

  Node xSquare = in[0]*in[0];
  Node ySquare = in[1]*in[1];
  Node zSquare = in[2]*in[2];

  Node s = sqrt((xSquare + ySquare + zSquare)*(xSquare + ySquare + zSquare)*(xSquare + ySquare + zSquare));
  Node trzy = Node(3.);
  Node jeden = Node(1.);

  out[3] = in[0]*(trzy*cs*(in[5]*in[1]-in[4]*in[2])-ce)/s;
  out[4] = (in[5]*(cs-trzy*cs*xSquare)    - ce*in[1] + trzy*cs*in[3]*in[0]*in[2])/s;
  out[5] = (cs*in[4]*(trzy*xSquare-jeden) - ce*in[2] - trzy*cs*in[3]*in[0]*in[1])/s;

  out[0] = in[3]+Node(0.);
  out[1] = in[4]+Node(0.);
  out[2] = in[5]+Node(0.);
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
  int dim = 4;

  // funkcja getV2 jest identyczna jak getV tylko jest uproszczona - os z nie jest uwzgledniana
  Map f(dim==4? getV2 : getV,dim,dim,paramsNumber);

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
  cout << "domain = [" << solution.getLeftDomain() << "," << solution.getRightDomain() << "]\n";
  cout << "------" << endl;

  int N = 200;
  //for(int i=0;i<=N;++i)
    //show(i*solution.getRightDomain()/N,solution,f);

  for(double d = 0.0 ; d < 30.5; d = d + 0.01){
    cout << "t: " << d << solution(d) << endl;
  }
  cout << "------" << endl;
}
