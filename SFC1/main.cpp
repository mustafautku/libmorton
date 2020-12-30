/*
 * main.cpp
 *
 *  Created on: Dec 28, 2020
 *      Author: mustafautku
 */




/*
 * main.cpp
 *
 *  Created on: Dec 15, 2020
 *      Author: mustafautku
 */
// Utility headers
#include "morton.h"
#include "math.h"
#include <iostream>
#include <spatialindex/SpatialIndex.h>

using namespace std;
using namespace libmorton;
using namespace SpatialIndex;

#define DEBUG 0
#define TESTING 0

static int DIMSCALE=16;  // #bits to represent single dim.
//std::map<uint32_t, uint32_t> resultset;
std::vector<uint32_t> resultset;
//uint32_t currmin;
//uint32_t currmax;

void findMortonRanges(Region *, Region *, uint32_t, int );

int main(int argc, char *argv[]) {

	// interleaving order (x,y)  : shape of z: z uclari yukarı aşagi bakıyor
	// interleaving order (y,x)  : shape of z: z gibi.

	// bu arada y,x sırasıyla interleave ediyor.

	Tools::Random rnd;
	int scale = pow(2, 16) - 1;
	if(TESTING){
		cout << morton2D_64_encode(0, 0) << endl;
		cout << morton2D_64_encode(1, 0) << endl;
		cout << morton2D_64_encode(50, 12) << endl;

		double x, y;

		x = 0.5;
		y = 0.5;
		int xs = x * scale;
		int ys = y * scale;
		cout << morton2D_32_encode(xs, ys) << endl;


		x = rnd.nextUniformDouble();
		y = rnd.nextUniformDouble();
		cout << x << ", " << y << endl;

		cout << morton2D_32_encode(x * scale, y * scale) << endl;

		// Test bit-wise operations
		uint32_t v1 = 1;
		cout << "v1: " << v1 << endl;	// 00001
		v1 = (v1 << 1);
		cout << "v1: " << v1 << endl;	// 00010
		v1 = (v1 << 1);
		cout << "v1: " << v1 << endl;	// 00100
		v1 = (v1 << 2);
		cout << "v1: " << v1 << endl;  // 10000

		uint32_t v2 = 3;
		cout << "v2: " << v2 << endl;  //11
		v2 = (v2 << 2);
		cout << "v2: " << v2 << endl;  // 1100

		uint32_t v3 = v1 & v2;  // 10000 & 1100 = 0
		uint32_t v4 = v1 | v2;  // 10000 | 1100 = 28
		cout << "v3: " << v3 << endl;
		cout << "v4: " << v4 << endl;

		set<int> myset;   // sorts data by default. :(  Does not insert duplicate :)
		myset.insert(10);
		myset.insert(8);
		myset.insert(17);
		myset.insert(11);
		myset.insert(11);
		set<int>::iterator itr;
		cout << "\nThe set s1 is : \n";  // 8,10,11,17,
		for (itr = myset.begin();
				itr != myset.end(); ++itr)
		{
			cout << *itr << ",";
		}
		cout << endl;

		vector<int> g1;
		g1.push_back(10);
		g1.push_back(8);
		g1.push_back(17);
		g1.push_back(11);
		g1.push_back(11);
		cout << "Output of vector: ";   // 10, 8, 17, 11, 11,
		for (auto i = g1.begin(); i != g1.end(); ++i)
			cout << *i << ", ";
		cout << endl;
	}

	// FINDING "morton-ranges" of a range query



	// TEST QUERY: (0,0) (2*1/2^16, 1/2^16) Test the SMALLEST cells.
//	double xl = 2*1/pow(2,16);
//	double yl = 0.0;//2*1/pow(2,16);
//	double xh = 5*1/pow(2,16);
//	double yh = 2*1/pow(2,16);
//	double *plow= new double[2]{xl,yl};
//	double *phigh= new double[2]{xh,yh};

	// TEST QUERY: (0,0) (2*1/2^16, 1/2^16) Test the LARGEST cells.
//	double xl = 1-1/pow(2,16);
//	double yl = 1-3*1/pow(2,16);
//	double xh = 1.0;
//	double yh = 1.0;
//	double *plow= new double[2]{xl,yl};
//	double *phigh= new double[2]{xh,yh};

	// RANDOM QUERY:
	double xl,xh;// = xl+dx;
	double yl,yh;// = yl+dy;

	do{   // generate a random query within unit area
		xl = rnd.nextUniformDouble();
		yl = rnd.nextUniformDouble();
		double dx = rnd.nextUniformDouble(0.3, 0.4);
		double dy = rnd.nextUniformDouble(0.3, 0.4);
		xh=xl+dx;
		yh=yl+dy;
	}while(xh>1.0 || yh>1.0);
	cout<< std::setprecision(25)<< xl << " "<<std::setprecision(25)<< xh << " "<<std::setprecision(25)<<yl << " "<<std::setprecision(25)<<yh << " "<<endl;

	// Range QUERY corners are aligned to the nearest multiple of (1/pow(2,16)). This approx. increases the exec. time without loosing correctness.
	int yuvarlama=(int)(xl/(1/pow(2,16)));
	xl=(1/pow(2,16))*yuvarlama;

	yuvarlama=(int)(xh/(1/pow(2,16)));
	xh=(1/pow(2,16))*yuvarlama;

	yuvarlama=(int)(yl/(1/pow(2,16)));
	yl=(1/pow(2,16))*yuvarlama;

	yuvarlama=(int)(yh/(1/pow(2,16)));
	yh=(1/pow(2,16))*yuvarlama;
	cout<< std::setprecision(25)<<xl << " "<<std::setprecision(25)<<xh << " "<<std::setprecision(25)<<yl << " "<<std::setprecision(25)<<yh << " "<<endl;

	double *plow= new double[2]{xl,yl};
	double *phigh= new double[2]{xh,yh};
	cout << morton2D_32_encode((uint16_t)(xl*scale),(uint16_t)(yl*scale)) << " -- ";
	cout << morton2D_32_encode((uint16_t)(xh*scale),(uint16_t)(yh*scale)) << endl;


	Region *query = new  Region(plow,phigh,2);
	cout << "RANGE QUERY: "<< *query<<endl;
	delete plow;
	delete phigh;

	// UNIT REGION:
	plow= new double[2]{0.0,0.0};
	phigh= new double[2]{1.0,1.0};
	Region *unit= new Region(plow,phigh,2);
	delete plow;
	delete phigh;
	int pos=1;
	uint32_t mc=0;
//	currmin=0;
//	currmax=0;
	findMortonRanges(query,unit,mc,pos);

	cout << "Output of vector: "<<endl;   // 10, 8, 17, 11, 11,
	for (auto i = resultset.begin(); i != resultset.end(); ++i)
		cout << *i << endl;
	cout << endl;

	// Convert to discrete partitions:
//	std::map<uint32_t,uint32_t> resultmap;
//	uint32_t* currmin=resultset.begin();
//	uint32_t* currmax=resultset.begin();
//	for (auto i = resultset.begin(); i != resultset.end(); ++i){
//		if(*i > currmax+1)
//
//	}

	delete query;
	delete unit;
}

void findMortonRanges(Region *rq, Region *runit, uint32_t mc,int pos){


	if(DEBUG) cout <<*runit<<endl;

	uint32_t quad0=(mc<<2);
	uint32_t quad1=quad0+1;
	uint32_t quad2=quad0+2;
	uint32_t quad3=quad0+3;



	Point pcenter;
	runit->getCenter(pcenter);

	{  // QUAD-0: LEFT-LOW SUBREGION
		double *p = new double[2] { runit->m_pLow[0], runit->m_pLow[1] };
		Point pll(p, 2);
		Region *subunit = new Region(pll, pcenter);
		if (DEBUG)
			cout << (*subunit) << endl;
		if (rq->containsRegion(*subunit)) {
			int shift = (2*DIMSCALE - 2 * pos);
			uint32_t mortonmin = (quad0 << shift);
			uint32_t mask = pow(2, shift) - 1;
			uint32_t mortonmax = (mortonmin | mask);

			resultset.push_back(mortonmin);
			resultset.push_back(0);  // a signal to indicate all values..
			resultset.push_back(mortonmax);

			//cout << mortonmin << " " << mortonmax<< ", ";
		} else if (rq->getIntersectingArea(*subunit) > 0) {
			if (pos < DIMSCALE)
				findMortonRanges(rq, subunit, quad0, pos + 1);
			else {
				if(DEBUG){
					int shift = (2*DIMSCALE - 2 * pos);
					uint32_t mortonmin = (quad0 << shift);
					uint32_t mask = pow(2, shift) - 1;
					uint32_t mortonmax = (mortonmin | mask);
					assert(mortonmin==mortonmax);
				}
				resultset.push_back(quad0);
				//cout << mortonmin << " " << mortonmax<<", ";
			}
		}

		delete p;
		delete subunit;
	}
	{  // QUAD-1: RIGHT-LOW SUBREGION
		double *p = new double[2] { (runit->m_pLow[0] + runit->m_pHigh[0]) / 2,
				runit->m_pLow[1] };
		Point pll(p, 2);
		double *pp = new double[2] { runit->m_pHigh[0], (runit->m_pLow[1]
				+ runit->m_pHigh[1]) / 2 };
		Point prh(pp, 2);
		Region *subunit = new Region(pll, prh);
		if (DEBUG)
			cout << *subunit << endl;
		if (rq->containsRegion(*subunit)) {
			int shift = (2*DIMSCALE  - 2 * pos);
			uint32_t mortonmin = (quad1 << shift);
			uint32_t mask = pow(2, shift) - 1;
			uint32_t mortonmax = (mortonmin | mask);
			resultset.push_back(mortonmin);
			resultset.push_back(0);  // a signal to indicate all values..
			resultset.push_back(mortonmax);
//			cout << mortonmin << " " << mortonmax << ", ";
		} else if (rq->getIntersectingArea(*subunit) > 0) {
			if (pos < DIMSCALE)
				findMortonRanges(rq, subunit, quad1, pos + 1);
			else {
				if(DEBUG){
					int shift = (2*DIMSCALE - 2 * pos);
					uint32_t mortonmin = (quad1 << shift);
					uint32_t mask = pow(2, shift) - 1;
					uint32_t mortonmax = (mortonmin | mask);
					assert(mortonmin==mortonmax);
				}
				resultset.push_back(quad1);
				//cout << mortonmin << " " << mortonmax<<", ";
			}
		}
		delete p;
		delete pp;
		delete subunit;
	}
	{  // QUAD-2: LEFT-HIGH SUBREGION
		double *p = new double[2] { runit->m_pLow[0], (runit->m_pLow[1]
				+ runit->m_pHigh[1]) / 2 };
		Point pll(p, 2);  //left-low-point
		double *pp = new double[2] { (runit->m_pLow[0] + runit->m_pHigh[0]) / 2,
				runit->m_pHigh[1] };
		Point prh(pp, 2);  //right-high-point
		Region *subunit = new Region(pll, prh);
		if (DEBUG)
			cout << *subunit << endl;
		if (rq->containsRegion(*subunit)) {
			int shift = (2*DIMSCALE  - 2 * pos);
			uint32_t mortonmin = (quad2 << shift);
			uint32_t mask = pow(2, shift) - 1;
			uint32_t mortonmax = (mortonmin | mask);
			resultset.push_back(mortonmin);
			resultset.push_back(0);  // a signal to indicate all values..
			resultset.push_back(mortonmax);
//			cout << mortonmin << " " << mortonmax << ", ";
		} else if (rq->getIntersectingArea(*subunit) > 0) {
			if (pos < DIMSCALE)
				findMortonRanges(rq, subunit, quad2, pos + 1);
			else {
				if(DEBUG){
					int shift = (2*DIMSCALE - 2 * pos);
					uint32_t mortonmin = (quad2 << shift);
					uint32_t mask = pow(2, shift) - 1;
					uint32_t mortonmax = (mortonmin | mask);
					assert(mortonmin==mortonmax);
				}
				resultset.push_back(quad2);
				//cout << mortonmin << " " << mortonmax<<", ";
			}
		}

		delete p;
		delete pp;
		delete subunit;
	}
	{  // QUAD-3: RIGTH-HIGH SUBREGION
		double *pp = new double[2] { runit->m_pHigh[0], runit->m_pHigh[1] };
		Point prh(pp, 2);  //right-high-point
		Region *subunit = new Region(pcenter, prh);
		if (DEBUG)
			cout << *subunit << endl;
		if (rq->containsRegion(*subunit)) {
			int shift = (2*DIMSCALE  - 2 * pos);
			uint32_t mortonmin = (quad3 << shift);
			uint32_t mask = pow(2, shift) - 1;
			uint32_t mortonmax = (mortonmin | mask);
			resultset.push_back(mortonmin);
			resultset.push_back(0);  // a signal to indicate all values..
			resultset.push_back(mortonmax);
//				cout << mortonmin << " " << mortonmax << ", ";
		} else if (rq->getIntersectingArea(*subunit) > 0) {
			if (pos < DIMSCALE)
				findMortonRanges(rq, subunit, quad3, pos + 1);
			else {
				if(DEBUG){
					int shift = (2*DIMSCALE - 2 * pos);
					uint32_t mortonmin = (quad3 << shift);
					uint32_t mask = pow(2, shift) - 1;
					uint32_t mortonmax = (mortonmin | mask);
					assert(mortonmin==mortonmax);
				}
				resultset.push_back(quad3);
				//cout << mortonmin << " " << mortonmax<<", ";
			}
		}

		delete pp;
		delete subunit;
	}
}
