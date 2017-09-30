
 #include <iostream>
#include <stdexcept>
#include "capd/capdlib.h"
#include "capd/krak/krak.h"
#include "capd/mpcapdlib.h"
#include <sstream>
#include <fstream>
#include "capd/poincare/NonlinearSection.h"
#include "capd/poincare/PoincareMap.h"
#include <cmath>

using namespace std;
using namespace capd;
using capd::autodiff::Node;




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

    sin(in[0]);

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

void find_and_replace(string& source, string const& find, string const& replace)
{
    for(string::size_type i = 0; (i = source.find(find, i)) != string::npos;)
    {
        source.replace(i, find.length(), replace);
        i += replace.length();
    }
}

string convert(long double myLongDouble) {
    stringstream blah;
    blah << myLongDouble;

    return blah.str();
}

void intro();
void part(string, string, string);
void partLast(string, string, string);
void ending();
void matlab_title(double, double, double, double);

int main(int argc, char* argv[]){

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

    int order = 10;
    OdeSolver s(f,order, OdeSolver::StepControlType(2,1e-12));
    s.setAbsoluteTolerance(1e-20);
    s.setRelativeTolerance(1e-20);

    Real initTime = 0.0;
    Real finalTime = 30;
    Real startTime = 0.0;

    TimeMap::SolutionCurve solution(initTime);
    TimeMap tm3(s);

    //istringstream ss1(argv[1]);
    //istringstream ss2(argv[2]);
    //istringstream ss3(argv[3]);
    istringstream ss4(argv[4]);

    double real1;
    double real2;
    double real3;
    double counter;





    //ss1 >> real1;
    //ss2 >> real2;
    //ss3 >> real3;
    ss4 >> counter;

    // Define Poincare section
    //LDNonlinearSection section("var:x,y,z;fun:z-1");
    //LDPoincareMap pm(s, section);

    istringstream ssR(argv[1]);
    istringstream ssPhi(argv[2]);
    istringstream ssTheta(argv[3]);


    double rValue = 1.0;
    double phiValue;
    double thetaValue = 2.61799387799;

    //ssR >> rValue;
    ssPhi >> phiValue;
    //ssTheta >> thetaValue;




    Real x(rValue * cos(thetaValue) * cos(phiValue));
    Real y(rValue * cos(thetaValue) * sin(phiValue));
    Real z(rValue * sin(thetaValue));


  Real d[] = {Real(x),Real(y),Real(z),Real(0.0),Real(.0),Real(.0)};
  //Real d[] = {Real(1.0),Real(1.0),Real(1.0),Real(0.0),Real(-2e-2),Real(2e-2)};


  Vector u3(dim,d);

  tm3(finalTime, u3,solution);

  int N = 500;

  //double finalTime = 300;
  double step = 0.01;

  double treshold = 29.99;



    ofstream myfile;
    //myfile.open("/home/marin/Desktop/project-capd/example2.txt");
    myfile.open(argv[5]);


    //matlab_title(real1, real2, real3, counter);
    myfile << "counter: " << counter <<  ", r: " << rValue << ", phi: " << phiValue << ", theta: " << thetaValue << "\n";




    for(double d = 0.0 ; d <= 30; d = d + step){
        //cout << "t: " << d << " : "<< solution(d) << endl;

        // plot 3d output
        string text0 = convert(solution(d)[0]);
        string text1 = convert(solution(d)[1]);
        string text2 = convert(solution(d)[2]);

        //find_and_replace(text0, "e", "*^");
        //find_and_replace(text1, "e", "*^");
        //find_and_replace(text2, "e", "*^");

        //cout << text0 << "," << text1 << "," << text2 << endl;
        myfile << text0 << "," << text1 << "," << text2 << "\n";


/*
/*
        if(d > treshold){
                cout << "{" << text0 << "," << text1 << "," << text2 << "}";
              //partLast(text0, text1, text2);



        } else {
               cout << "{" << text0 << "," << text1 << "," << text2 << "}, ";
             //part(text0, text1, text2);

        }
        */


        //cout << "t: " << d << " : "<< solution(d)[0] << " , "<< solution(d)[1] << " , "<< solution(d)[2] << endl;
    }

    myfile.close();


}

void matlab_title(double d1, double d2, double d3, double counter){
    cout << "counter: " << counter <<  ", vx: " << d1 << ", vy: " << d2 << ", vz: " << d3 << "\n";

}

void intro(){
    cout << "(* Content-type: application/vnd.wolfram.mathematica *)" << endl;
    cout << ""<< endl;
    cout << "(*** Wolfram Notebook File ***)"<< endl;
    cout << "(* http://www.wolfram.com/nb *)"<< endl;
    cout << ""<< endl;
    cout << "(* CreatedBy='Mathematica 10.3' *)"<< endl;
    cout << ""<< endl;
    cout << "(*CacheID: 234*)"<< endl;
    cout << "(* Internal cache information:"<< endl;
    cout << "NotebookFileLineBreakTest"<< endl;
    cout << "NotebookFileLineBreakTest"<< endl;
    cout << "NotebookDataPosition[       158,          7]"<< endl;
    cout << "NotebookDataLength[      1390,         50]"<< endl;
    cout << "NotebookOptionsPosition[      1151,         37]"<< endl;
    cout << "NotebookOutlinePosition[      1494,         52]"<< endl;
    cout << "CellTagsIndexPosition[      1451,         49]"<< endl;
    cout << "WindowFrame->Normal*)"<< endl;
    cout << ""<< endl;
    cout << "(* Beginning of Notebook Content *)"<< endl;
    cout << "Notebook[{"<< endl;
    cout << "Cell[BoxData["<< endl;
    cout << "RowBox[{\"data\", \"=\", "<< endl;
    cout << "RowBox[{\"{\","<< endl;
    cout << "RowBox[{"<< endl;
}

void part(string d1, string d2, string d3){
    cout << "RowBox[{\"{\"," << endl;
    cout << "RowBox[{\"" + d1 + "\", \",\", \"" + d2 + "\", \",\", \"" + d3 + "\"}], \"}\"}], \",\", "<< endl;
}

void partLast(string d0, string d1, string d2){
    cout << "RowBox[{\"{\"," << endl;

    cout << "RowBox[{\"" + d0 + "\", \",\", \"" + d1 + "\", \",\", \"" + d2 + "\"}], \"}\"}]}], " << endl;
}

void ending(){
    cout << "\"}\"}]}]], \"Input\"," << endl;
    cout << "CellChangeTimes->{{3.7229380911609364`*^9, 3.7129380928383675`*^9}}]" << endl;
    cout << "},"<< endl;
    cout << "WindowSize->{1510, 781}," << endl;
    cout << "WindowMargins->{{0, Automatic}, {Automatic, 0}}," << endl;
    cout << "FrontEndVersion->\"10.3 for Microsoft Windows (64-bit) (October 9, 2015)\","<< endl;
    cout << "StyleDefinitions->\"Default.nb\""<< endl;
    cout << "]"<< endl;
    cout << "(* End of Notebook Content *)"<< endl;
    cout << ""<< endl;
    cout << "(* Internal cache information *)"<< endl;
    cout << "(*CellTagsOutline"<< endl;
    cout << "CellTagsIndex->{}"<< endl;
    cout << "*)"<< endl;
    cout << "(*CellTagsIndex"<< endl;
    cout << "CellTagsIndex->{}"<< endl;
    cout << "*)"<< endl;
    cout << "(*NotebookFileOutline"<< endl;
    cout << "Notebook[{"<< endl;
    cout << "Cell[558, 20, 589, 15, 31, \"Input\"]"<< endl;
    cout << "}"<< endl;
    cout << "]"<< endl;
    cout << "*)"<< endl;
    cout << ""<< endl;
    cout << "(* End of internal cache information *)"<< endl;
    cout << ""<< endl;
}


