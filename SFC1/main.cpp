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

static int DIMSCALE=15;
std::map<uint32_t, uint32_t> resultset;

void findMortonRanges(Region *, Region *, uint32_t, int );

int main(int argc, char *argv[]) {

	// interleaving order (x,y)  : shape of z: z uclari yukarı aşagi bakıyor
	// interleaving order (y,x)  : shape of z: z gibi.

	// bu arada y,x sırasıyla interleave ediyor.

	Tools::Random rnd;
	if(TESTING){
		cout << morton2D_64_encode(0, 0) << endl;
		cout << morton2D_64_encode(1, 0) << endl;
		cout << morton2D_64_encode(50, 12) << endl;

		double x, y;
		int scale = pow(2, 16) - 1;
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
	}
	// FINDING "morton-ranges" of a range query
//
	double x = rnd.nextUniformDouble();
	double y = rnd.nextUniformDouble();
	double dx = rnd.nextUniformDouble(0.01, 0.1);
	double dy = rnd.nextUniformDouble(0.01, 0.1);
	double *plow= new double[2]{x,y};
	double *phigh= new double[2]{x+dx,y+dy};

	// QUERY:
//	double *plow= new double[2]{1-1/pow(2,16),0.0};
//	double *phigh= new double[2]{1,1/pow(2,16)};
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

	findMortonRanges(query,unit,mc,pos);

	delete query;
	delete unit;
}

void findMortonRanges(Region *rq, Region *runit, uint32_t mc,int pos){


	if(DEBUG) cout <<*runit<<endl;

	uint32_t quad0=(mc<<2);
	uint32_t quad1=quad0+1;
	uint32_t quad2=quad0+2;
	uint32_t quad3=quad0+3;

	uint32_t currmin;
	uint32_t currmax;

	Point pcenter;
	runit->getCenter(pcenter);

	{  // QUAD-0: LEFT-LOW SUBREGION
		double *p=new double[2]{runit->m_pLow[0],runit->m_pLow[1]};
		Point pll(p,2);
		Region *subunit= new Region(pll,pcenter);
		if(DEBUG) cout<<(*subunit)<<endl;
		if(rq->containsRegion(*subunit)){
			int shift= (32-2*pos);
			uint32_t mortonmin = (quad0<<shift);
			uint32_t mask=pow(2,shift)-1;
			uint32_t mortonmax = (mortonmin|mask);
			resultset.insert(std::pair<uint32_t, uint32_t>(mortonmin, mortonmax));
			currmin=mortonmin;
			currmax=mortonmax;
			//cout << mortonmin << " " << mortonmax<< ", ";
		}
		else if(rq->getIntersectingArea(*subunit)>0){
			if(pos<DIMSCALE)
				findMortonRanges(rq,subunit,quad0,pos+1);
			else{
				int shift= (32-2*pos);
				uint32_t mortonmin = (quad0<<shift);
				uint32_t mask=pow(2,shift)-1;
				uint32_t mortonmax = (mortonmin|mask);
				resultset.insert(std::pair<uint32_t, uint32_t>(mortonmin, mortonmax));
				currmin=mortonmin;
				currmax=mortonmax;
				//cout << mortonmin << " " << mortonmax<<", ";
			}
		}

		delete p;
		delete subunit;
	}
	{// QUAD-1: RIGHT-LOW SUBREGION
		double *p = new double[2] { (runit->m_pLow[0] +runit->m_pHigh[0])/2, runit->m_pLow[1] };
		Point pll(p, 2);
		double *pp = new double[2] { runit->m_pHigh[0], (runit->m_pLow[1]+ runit->m_pHigh[1])/2 };
		Point prh(pp, 2);
		Region *subunit = new Region(pll, prh);
		if(DEBUG) cout<<*subunit<<endl;
		if (rq->containsRegion(*subunit)) {
			int shift = (32 - 2 * pos);
			uint32_t mortonmin = (quad1 << shift);
			uint32_t mask = pow(2, shift) - 1;
			uint32_t mortonmax = (mortonmin | mask);
			if(mortonmin==currmax+1){
				currmax=mortonmax;
				resultset.erase(currmin);
				resultset.insert(std::pair<uint32_t, uint32_t>(currmin, currmax));
			}
			else{
				resultset.insert(std::pair<uint32_t, uint32_t>(currmin, currmax));
				currmin=mortonmin;
				currmax=mortonmax;
			}
//			cout << mortonmin << " " << mortonmax << ", ";
		} else if(rq->getIntersectingArea(*subunit)>0){
			if(pos<DIMSCALE)
				findMortonRanges(rq, subunit, quad1, pos+1);
			else{
				int shift= (32-2*pos);
				uint32_t mortonmin = (quad1<<shift);
				uint32_t mask=pow(2,shift)-1;
				uint32_t mortonmax = (mortonmin|mask);
				if(mortonmin==currmax+1)
					currmax=mortonmax;
				else{
					resultset.insert(std::pair<uint32_t, uint32_t>(currmin, currmax));
					currmin=mortonmin;
					currmax=mortonmax;
				}
//				cout << mortonmin << " " << mortonmax<<", ";
			}
		}
		delete p;
		delete pp;
		delete subunit;
	}
	{// QUAD-2: LEFT-HIGH SUBREGION
		double *p = new double[2] { runit->m_pLow[0], (runit->m_pLow[1]+runit->m_pHigh[1])/2 };
		Point pll(p, 2);  //left-low-point
		double *pp = new double[2] { (runit->m_pLow[0]+runit->m_pHigh[0])/2, runit->m_pHigh[1]};
		Point prh(pp, 2);  //right-high-point
		Region *subunit = new Region(pll, prh);
		if(DEBUG) cout<<*subunit<<endl;
		if (rq->containsRegion(*subunit)) {
			int shift = (32 - 2 * pos);
			uint32_t mortonmin = (quad2 << shift);
			uint32_t mask = pow(2, shift) - 1;
			uint32_t mortonmax = (mortonmin | mask);
			if(mortonmin==currmax+1)
				currmax=mortonmax;
			else{
				resultset.insert(std::pair<uint32_t, uint32_t>(currmin, currmax));
				currmin=mortonmin;
				currmax=mortonmax;
			}
//			cout << mortonmin << " " << mortonmax << ", ";
		} else if (rq->getIntersectingArea(*subunit) > 0) {
			if (pos<DIMSCALE)
				findMortonRanges(rq, subunit, quad2, pos + 1);
			else {
				int shift = (32 - 2 * pos);
				uint32_t mortonmin = (quad2 << shift);
				uint32_t mask = pow(2, shift) - 1;
				uint32_t mortonmax = (mortonmin | mask);
				if(mortonmin==currmax+1)
					currmax=mortonmax;
				else{
					resultset.insert(std::pair<uint32_t, uint32_t>(currmin, currmax));
					currmin=mortonmin;
					currmax=mortonmax;
				}
//				cout << mortonmin << " " << mortonmax << ", ";
			}
		}

		delete p;
		delete pp;
		delete subunit;
	}
	{// QUAD-3: RIGTH-HIGH SUBREGION
			double *pp = new double[2] { runit->m_pHigh[0], runit->m_pHigh[1]};
			Point prh(pp, 2);  //right-high-point
			Region *subunit = new Region(pcenter, prh);
			if(DEBUG) cout<<*subunit<<endl;
			if (rq->containsRegion(*subunit)) {
				int shift = (32 - 2 * pos);
				uint32_t mortonmin = (quad3 << shift);
				uint32_t mask = pow(2, shift) - 1;
				uint32_t mortonmax = (mortonmin | mask);
				if(mortonmin==currmax+1)
					currmax=mortonmax;
				else{
					resultset.insert(std::pair<uint32_t, uint32_t>(currmin, currmax));
					currmin=mortonmin;
					currmax=mortonmax;
				}
//				cout << mortonmin << " " << mortonmax << ", ";
			} else if(rq->getIntersectingArea(*subunit)>0){
				if(pos<DIMSCALE)
					findMortonRanges(rq, subunit, quad3, pos+1);
				else{
					int shift= (32-2*pos);
					uint32_t mortonmin = (quad3<<shift);
					uint32_t mask=pow(2,shift)-1;
					uint32_t mortonmax = (mortonmin|mask);
					if(mortonmin==currmax+1)
						currmax=mortonmax;
					else{
						resultset.insert(std::pair<uint32_t, uint32_t>(currmin, currmax));
						currmin=mortonmin;
						currmax=mortonmax;
					}
//					cout << mortonmin << " " << mortonmax<<", ";
				}
			}

			delete pp;
			delete subunit;
		}
}
