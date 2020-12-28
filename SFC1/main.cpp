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


//class Region
//{
//public:
//	double m_xmin, m_ymin, m_xmax, m_ymax;
//
//	Region(double x1, double y1, double x2, double y2)
//	{
//		m_xmin = (x1 < x2) ? x1 : x2;
//		m_ymin = (y1 < y2) ? y1 : y2;
//		m_xmax = (x1 > x2) ? x1 : x2;
//		m_ymax = (y1 > y2) ? y1 : y2;
//	}
//};

void findMortonRanges(Region *, Region *, uint32_t, int );

int main(int argc, char *argv[]) {

	// interleaving order (x,y)  : shape of z: z uclari yukarı aşagi bakıyor
	// interleaving order (y,x)  : shape of z: z gibi.

	// bu rada y,x sırasıyla interleave ediyor.
	cout<<morton2D_64_encode(0,0) << endl;

	cout<<morton2D_64_encode(1,0) << endl;
	cout<<morton2D_64_encode(50,12) << endl;

	double x,y;
	int scale=pow(2,16)-1;
	x=0.5;
	y=0.5;
	int xs=x*scale;
	int ys=y*scale;
	cout<<morton2D_32_encode(xs,ys) <<endl;

	Tools::Random rnd;
	x=rnd.nextUniformDouble();
	y=rnd.nextUniformDouble();
	cout << x <<", " << y  <<endl;

	cout<<morton2D_32_encode(x*scale,y*scale) <<endl;

	// Test bit-wise operations
	uint32_t v1=1;
	cout <<"v1: "<< v1<< endl;// 00001
	v1= (v1<<1);
	cout <<"v1: "<< v1<< endl;// 00010
	v1= (v1<<1);
	cout <<"v1: "<< v1<< endl;// 00100
	v1= (v1<<2);
	cout <<"v1: "<< v1<< endl;  // 10000

	uint32_t v2=3;
	cout <<"v2: "<< v2<< endl;  //11
	v2= (v2<<2);
	cout <<"v2: "<< v2<< endl;  // 1100

	uint32_t v3= v1&v2;  // 10000 & 1100 = 0
	uint32_t v4= v1|v2;  // 10000 | 1100 = 28
	cout <<"v3: "<< v3<< endl;
	cout <<"v4: "<< v4<< endl;

	// FINDING "morton-ranges" of a range query
//
//	x = rnd.nextUniformDouble();
//	y = rnd.nextUniformDouble();
//	double dx = rnd.nextUniformDouble(0.01, 0.1);
//	double dy = rnd.nextUniformDouble(0.01, 0.1);
//	double *plow= new double[2]{x,y};
//	double *phigh= new double[2]{x+dx,y+dy};
	double *plow= new double[2]{0.25,0.0};
	double *phigh= new double[2]{0.75,0.25};
	Region *q = new  Region(plow,phigh,2);
	delete plow;
	delete phigh;

	plow= new double[2]{0.0,0.0};
	phigh= new double[2]{1.0,1.0};
	Region *unit= new Region(plow,phigh,2);
	delete plow;
	delete phigh;
	int pos=1;
	uint32_t mc=0;
	findMortonRanges(q,unit,mc,pos);
	delete unit;
}

void findMortonRanges(Region *rq, Region *runit, uint32_t mc,int pos){


	cout <<*runit<<endl;

	uint32_t quad0=(mc<<2);
	uint32_t quad1=quad0+1;
	uint32_t quad2=quad0+2;
	uint32_t quad3=quad0+3;

	Point p11;
	runit->getCenter(p11);

	double *p=new double[2]{runit->m_pLow[0],runit->m_pLow[1]};
	Point p00(p,2);
	Region *runit0= new Region(p00,p11);
	cout<<(*runit0)<<endl;
	if(rq->containsRegion(*runit0)){
		int shift= (32-2*pos);
		uint32_t quad0min = (quad0<<shift);
		uint32_t mask=pow(2,shift)-1;
		uint32_t quad0max = (quad0min|mask);
		cout << quad0min << " " << quad0max<<endl;
	}
	else if(rq->getIntersectingArea(*runit0)>0)
		findMortonRanges(rq,runit0,quad0,pos+1);

	delete p;

	p = new double[2] { runit->m_pHigh[0]/2, runit->m_pLow[1] };
	Point p01(p, 2);
	double *pp = new double[2] { runit->m_pHigh[0], runit->m_pHigh[1]/2 };
	Point p12(pp, 2);
	Region *runit1 = new Region(p01, p12);
	cout<<*runit1<<endl;
	if (rq->containsRegion(*runit1)) {
		int shift = (32 - 2 * pos);
		uint32_t quad1min = (quad1 << shift);
		uint32_t mask = pow(2, shift) - 1;
		uint32_t quad1max = (quad1min | mask);
		cout << quad1min << " " << quad1max << endl;
	} else if(rq->getIntersectingArea(*runit1)>0)
		findMortonRanges(rq, runit1, quad1, pos+1);


}
