/*
 * main.cpp
 *
 *  Generated on: Dec 28, 2020
 *      Author: mustafautku
 *
 *
 *     Program finds RANGE QUERY's morton codes for different scalings, at most 32-bit mortons.
 *     I am using Quadtree technique to genereate the codes.
 *     to generate morton partitions of the query region, 16-bit- morton codes seems enough. That means DIMSCALE to be 8.
 *     For higher DIMSCALE values, the number of parititions increase very much.
 *     If you are going to apply the query (that is represented 16-bit morton partitions) on a a B-tree that is storing 32-bit spatial data morton codes,
 *     then you should convert the generated partitions to 32-bit-scale of mortons.
 *
 */

/*
 * main.cpp
 *
 *  Generated on: Dec 15, 2020
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
#define TESTING 1

static int QUERYDIMSCALE=8;  // DEPTH of z-curve = NUMBER OF MAX SPLITS = bits used to represent single dim.  8 seems enough to represent the shape of query.
// naming is a little confusing. But, the real meaning is we have a double in unit area [0,1]. If we want to represent the space 2^16 = 65535 cells (codes),
// each dim will have 2^8=256 slices. So, we have to *scale* double value to this range [0-2^8]
// Boylece unit area yı 65535 hucre ile temsil ettik.
// Precision daha iyi olmasini istersek DEPTH= MAX_SPLIT= 32 yapabilriz. Yani unit area yı 2^64 tane hucre ile temsile edeceğiz...

//static int MAXNUMCELLINDIM = pow(2, QUERYDIMSCALE) - 1;
std::map<uint32_t, uint32_t> resultset;
uint32_t currmin;
uint32_t currmax;

void findMortonRanges(Region *, Region *, uint32_t, int );

int main(int argc, char *argv[]) {

	// interleaving order (x,y)  : shape of z: z uclari yukarı aşagi bakıyor
	// interleaving order (y,x)  : shape of z: z gibi.

	// libmorton y,x sırasıyla interleave ediyor.

	Tools::Random rnd;


	if(TESTING){
		// Her eksen [0- 2^32] aralıgında bir deger bekliyor.
		cout << morton2D_64_encode(0, 0) << endl;
		cout << morton2D_64_encode(1, 0) << endl;  // y:00 x:01  --> 0001
		cout << morton2D_64_encode(50, 12) << endl; // y:12 x:50  --> 001100  110010  --> 010110100100 = 1444

		double x, y;

		x = 0.5;
		y = 0.5;
		int xs = x * (pow(2, QUERYDIMSCALE) - 1); // scale to [0,255]
		int ys = y * (pow(2, QUERYDIMSCALE) - 1); // scale to [0,255]

		cout << "max morcod:"<< morton2D_32_encode(255, 255)<< endl;  // 65535
		cout << "max morcod:"<< morton2D_64_encode(255, 255)<< endl;  // 65535. evet aynı deger. 65535 tane kod kullanmak istedik. Her iki calsma uzayı da bize bunu verdi.

		uint32_t morcod = morton2D_32_encode(xs, ys);
		cout << "morcod of center:"<< morcod<< endl;

		cout << "morcod:"<< morton2D_32_encode(256, 256)<< endl;  // 65535*3 =196608
		cout << "morcod:"<< morton2D_64_encode(256, 256)<< endl;  // 65535*3 =196608. evet aynı deger. 65535 tane kod kullanmak istedik. Her iki calsma uzayı da bize bunu verdi.

		// Yani uretilen morton kodlari sol alt kose pix'in morton kodu. Cell ne kadar buyuk olursa olsun aynı cunku.!!!

		uint_fast16_t x1,y1;
		morton2D_32_decode(morcod,x1,y1);
		cout << "decoding:" << x1 << '\t' << y1 << endl;  // shows the number: 2^16  / 2 = 32767
		cout << "decoding:" << (double)  x1/pow(2, QUERYDIMSCALE) << '\t' <<(double) y1/pow(2, QUERYDIMSCALE) << endl;  // move it into unit area.


		x = rnd.nextUniformDouble();
		y = rnd.nextUniformDouble();
		cout << x << ", " << y << endl;

		cout << morton2D_32_encode(x * (pow(2, QUERYDIMSCALE) - 1), y * (pow(2, QUERYDIMSCALE) - 1)) << endl;

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
		return 0;
	}



	// FINDING "morton-ranges" of a range query by using QuadTree approach.

	// TEST QUERY: (0,0) (2*1/2^16, 1/2^16) Test the SMALLEST cells.
//	double xl = 0.0; // 2*1/pow(2,16);
//	double yl = 0.0; // 2*1/pow(2,16);
//	double xh = 1.0; // 5*1/pow(2,16);
//	double yh = 1.0; //2*1/pow(2,16);
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
	cout<< std::setprecision(25) << xl << " "<<std::setprecision(25)<< xh << " "<<std::setprecision(25)<< yl << " "<<std::setprecision(25)<< yh << " "<<endl;

	// ALIGNMENT of query corners. Approximate doubles to multiples of 1/pow(2,16). this increase the exec. time without lose of accuracy.
	// thanks to this alignment in the findMortonRanges, it never enters "else if (rq->getIntersectingArea(*subunit) > 0) { else pos==QUERYDIMSCALE" section..
	int yuvarlama=(int)(xl/(1/pow(2,QUERYDIMSCALE)));
	xl=(1/pow(2,QUERYDIMSCALE))*yuvarlama;

	yuvarlama=(int)(xh/(1/pow(2,QUERYDIMSCALE)));
	xh=(1/pow(2,QUERYDIMSCALE))*yuvarlama;

	yuvarlama=(int)(yl/(1/pow(2,QUERYDIMSCALE)));
	yl=(1/pow(2,QUERYDIMSCALE))*yuvarlama;

	yuvarlama=(int)(yh/(1/pow(2,QUERYDIMSCALE)));
	yh=(1/pow(2,QUERYDIMSCALE))*yuvarlama;

	cout<<std::setprecision(25)<<  xl << " "<<std::setprecision(25)<< xh << " "<<std::setprecision(25)<< yl << " "<<std::setprecision(25)<< yh << " "<<endl;

	// Test case:  All values below are multiples of 1/pow(2,16), for example.
//	xl=0.0445404052734375;
//	xh=0.3786163330078125;
//	yl=0.100006103515625;
//	yh=0.45556640625;

	double *plow = new double[2]{xl,yl};
	double *phigh= new double[2]{xh,yh};


	cout << morton2D_32_encode((uint16_t)(xl*(pow(2, QUERYDIMSCALE) - 1)),(uint16_t)(yl*(pow(2, QUERYDIMSCALE) - 1))) << " -- ";
	cout << morton2D_32_encode((uint16_t)(xh*(pow(2, QUERYDIMSCALE) - 1)),(uint16_t)(yh*(pow(2, QUERYDIMSCALE) - 1))) << endl;


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

	int pos=1;  // start using 1 bit for each dim.
	uint32_t mc=0;
	currmin=0;
	currmax=0;

	findMortonRanges(query,unit,mc,pos);  // here traverse the quad-tree recursively

	// last partititon is being added below.
	if(resultset.find(currmin) != resultset.end()){
		resultset.erase(currmin);
		resultset.insert(
				std::pair<uint32_t, uint32_t>(currmin, currmax));
	}else
	resultset.insert(
					std::pair<uint32_t, uint32_t>(currmin, currmax));

	// scan the partitions
	cout <<"Number of partititions: "<< resultset.size()<<endl;
	for (auto itr = resultset.begin(); itr != resultset.end();itr++) {
		cout << '\t' << itr->first
				<< '\t' << itr->second << '\n';
	}
	cout << endl;


	cout << " POINTS in the QUERY RANGE:" << endl;
	// Now generate "double end points of partitions" to see the queryrange in gnuplot..PLOTS the range query correct. No error found.
	// generates many points at the edges and few points for the internal area. thus I made reduction by eliminating single "points(cells)"
	int reducedPartitions=0;
	for (auto itr = resultset.begin(); itr != resultset.end(); ++itr) {

		if( itr->second > itr->first  ){  //Remove smallest quads.  only select clustered cells(quads) or half of quad.
			reducedPartitions++;
			uint_fast16_t x1,y1, x2,y2;
			morton2D_32_decode(itr->first,x1,y1);
			morton2D_32_decode(itr->second,x2,y2);

			double xx1=(double)x1/pow(2, QUERYDIMSCALE);
			double yy1=(double)y1/pow(2, QUERYDIMSCALE);
			double xx2=(double)x2/pow(2, QUERYDIMSCALE);
			double yy2=(double)y2/pow(2, QUERYDIMSCALE);

			cout << std::setprecision(10) << xx1 << '\t' << std::setprecision(10)<< yy1 << endl;
			cout << std::setprecision(10) << xx2 << '\t' <<std::setprecision(10) << yy2 << endl;

			// convert 16-bit morton code to 32-bit codes:  2 ways: 1-transform by scaling over morton space 2- transform for each double dimension
//			uint32_t bit32mc1=( itr->first / pow(2, 2*QUERYDIMSCALE) ) * pow(2, 4* QUERYDIMSCALE);
//			uint32_t bit32mc2=( itr->second / pow(2,2*QUERYDIMSCALE) ) * pow(2, 4* QUERYDIMSCALE);
//
//			uint32_t bit32mc11= morton2D_32_encode(xx1*pow(2, 2*QUERYDIMSCALE), yy1*pow(2, 2*QUERYDIMSCALE));
//			uint32_t bit32mc22=morton2D_32_encode(xx2*pow(2, 2*QUERYDIMSCALE), yy2*pow(2, 2*QUERYDIMSCALE));
//
//			cout << bit32mc1 << ',' << bit32mc2 << endl;
//			cout << bit32mc11 << ',' << bit32mc22 << endl;

		}
	}
	cout <<"Number of partititions (after reduction): "<< reducedPartitions <<endl;  // Sometimes no reduction may not be possible depending on the query's locaiton.
		cout << endl;

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
			int shift = (2*QUERYDIMSCALE - 2 * pos);
			uint32_t mortonmin = (quad0 << shift);
			uint32_t mask = pow(2, shift) - 1;
			uint32_t mortonmax = (mortonmin | mask);

			if (resultset.empty()){
				currmin=mortonmin;
				currmax = mortonmax;
				resultset.insert(
						std::pair<uint32_t, uint32_t>(currmin, currmax));
			}
			else if (mortonmin == currmax + 1) {
				currmax = mortonmax;
			} else {

				if(resultset.find(currmin) != resultset.end()){
					resultset.erase(currmin);
					resultset.insert(
							std::pair<uint32_t, uint32_t>(currmin, currmax));
				}else
					resultset.insert(
							std::pair<uint32_t, uint32_t>(currmin, currmax));
				currmin = mortonmin;
				currmax = mortonmax;
			}
			//cout << mortonmin << " " << mortonmax<< ", ";
		} else if (rq->getIntersectingArea(*subunit) > 0) {
			if (pos < QUERYDIMSCALE)
				findMortonRanges(rq, subunit, quad0, pos + 1);
			else {// NEVER ENTERs HERE IFF CORNERS OF RANGE QUERY is ALIGNED TO MULTIPLES OF 1/pow(2,16). QUERYDIMSCALE=16
				int shift = (2*QUERYDIMSCALE - 2 * pos);
				uint32_t mortonmin = (quad0 << shift);
				uint32_t mask = pow(2, shift) - 1;
				uint32_t mortonmax = (mortonmin | mask);

				if (resultset.empty()){
					currmin=mortonmin;
					currmax = mortonmax;
					resultset.insert(
							std::pair<uint32_t, uint32_t>(currmin, currmax));
				}
				else if (mortonmin == currmax + 1) {
					currmax = mortonmax;
				} else {
					if(resultset.find(currmin) != resultset.end()){
						resultset.erase(currmin);
						resultset.insert(
								std::pair<uint32_t, uint32_t>(currmin, currmax));
					}else
						resultset.insert(
								std::pair<uint32_t, uint32_t>(currmin, currmax));
					currmin = mortonmin;
					currmax = mortonmax;
				}
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
			int shift = (2*QUERYDIMSCALE - 2 * pos);
			uint32_t mortonmin = (quad1 << shift);
			uint32_t mask = pow(2, shift) - 1;
			uint32_t mortonmax = (mortonmin | mask);
			if (resultset.empty()){
				currmin = mortonmin;
				currmax = mortonmax;
				resultset.insert(
										std::pair<uint32_t, uint32_t>(currmin, currmax));
			}else if (mortonmin == currmax + 1) {
				currmax = mortonmax;
			} else {
				resultset.insert(
						std::pair<uint32_t, uint32_t>(currmin, currmax));
				currmin = mortonmin;
				currmax = mortonmax;
			}
//			cout << mortonmin << " " << mortonmax << ", ";
		} else if (rq->getIntersectingArea(*subunit) > 0) {
			if (pos < QUERYDIMSCALE)
				findMortonRanges(rq, subunit, quad1, pos + 1);
			else {// NEVER ENTERs HERE IFF CORNERS OF RANGE QUERY is ALIGNED TO MULTIPLES OF 1/pow(2,16). QUERYDIMSCALE=16
				int shift = (2*QUERYDIMSCALE - 2 * pos);
				uint32_t mortonmin = (quad1 << shift);
				uint32_t mask = pow(2, shift) - 1;
				uint32_t mortonmax = (mortonmin | mask);
				if (resultset.empty()){
					currmin = mortonmin;
					currmax = mortonmax;
					resultset.insert(
							std::pair<uint32_t, uint32_t>(currmin, currmax));
				}else if (mortonmin == currmax + 1) {
					currmax = mortonmax;
				} else {
					resultset.insert(
							std::pair<uint32_t, uint32_t>(currmin, currmax));
					currmin = mortonmin;
					currmax = mortonmax;
				}
//				cout << mortonmin << " " << mortonmax<<", ";
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
			int shift = (2*QUERYDIMSCALE - 2 * pos);
			uint32_t mortonmin = (quad2 << shift);
			uint32_t mask = pow(2, shift) - 1;
			uint32_t mortonmax = (mortonmin | mask);
			if (resultset.empty()){
				currmin = mortonmin;
				currmax = mortonmax;
				resultset.insert(
						std::pair<uint32_t, uint32_t>(currmin, currmax));
			}else if (mortonmin == currmax + 1) {
				currmax = mortonmax;
			} else {
				if(resultset.find(currmin) != resultset.end()){
					resultset.erase(currmin);
					resultset.insert(
							std::pair<uint32_t, uint32_t>(currmin, currmax));
				}else
					resultset.insert(
							std::pair<uint32_t, uint32_t>(currmin, currmax));
				currmin = mortonmin;
				currmax = mortonmax;
			}
//			cout << mortonmin << " " << mortonmax << ", ";
		} else if (rq->getIntersectingArea(*subunit) > 0) {
			if (pos < QUERYDIMSCALE)
				findMortonRanges(rq, subunit, quad2, pos + 1);
			else { // NEVER ENTERs HERE IFF CORNERS OF RANGE QUERY is ALIGNED TO MULTIPLES OF 1/pow(2,16). QUERYDIMSCALE=16
				int shift = (2*QUERYDIMSCALE - 2 * pos);
				uint32_t mortonmin = (quad2 << shift);
				uint32_t mask = pow(2, shift) - 1;
				uint32_t mortonmax = (mortonmin | mask);
				if (resultset.empty()){
					currmin = mortonmin;
					currmax = mortonmax;
					resultset.insert(
							std::pair<uint32_t, uint32_t>(currmin, currmax));
				}else if(mortonmin == currmax + 1) {
					currmax = mortonmax;
				} else {
					if(resultset.find(currmin) != resultset.end()){
						resultset.erase(currmin);
						resultset.insert(
								std::pair<uint32_t, uint32_t>(currmin, currmax));
					}else
						resultset.insert(
								std::pair<uint32_t, uint32_t>(currmin, currmax));
					currmin = mortonmin;
					currmax = mortonmax;
				}
				//				cout << mortonmin << " " << mortonmax << ", ";
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
			int shift = (2*QUERYDIMSCALE - 2 * pos);
			uint32_t mortonmin = (quad3 << shift);
			uint32_t mask = pow(2, shift) - 1;
			uint32_t mortonmax = (mortonmin | mask);
			if (resultset.empty()){
				currmin = mortonmin;
				currmax = mortonmax;
				resultset.insert(
						std::pair<uint32_t, uint32_t>(currmin, currmax));
			}else if (mortonmin == currmax + 1) {
				currmax = mortonmax;
			} else {
				resultset.insert(
						std::pair<uint32_t, uint32_t>(currmin, currmax));
				currmin = mortonmin;
				currmax = mortonmax;
			}
//				cout << mortonmin << " " << mortonmax << ", ";
		} else if (rq->getIntersectingArea(*subunit) > 0) {
			if (pos < QUERYDIMSCALE)
				findMortonRanges(rq, subunit, quad3, pos + 1);
			else {  // NEVER ENTERs HERE IFF CORNERS OF RANGE QUERY is ALIGNED TO MULTIPLES OF 1/pow(2,16). QUERYDIMSCALE=16
				int shift = (2*QUERYDIMSCALE - 2 * pos);
				uint32_t mortonmin = (quad3 << shift);
				uint32_t mask = pow(2, shift) - 1;
				uint32_t mortonmax = (mortonmin | mask);
				if (resultset.empty()){
					currmin = mortonmin;
					currmax = mortonmax;
					resultset.insert(
							std::pair<uint32_t, uint32_t>(currmin, currmax));
				}else if (mortonmin == currmax + 1) {
					currmax = mortonmax;
				} else {
					resultset.insert(
							std::pair<uint32_t, uint32_t>(currmin, currmax));
					currmin = mortonmin;
					currmax = mortonmax;
				}
//					cout << mortonmin << " " << mortonmax<<", ";
			}
		}

		delete pp;
		delete subunit;
	}
}
