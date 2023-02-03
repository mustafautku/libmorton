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
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

using namespace std;
using namespace libmorton;
using namespace SpatialIndex;

#define DEBUG 0
#define TESTING 1

static int WIDTH=8;  // DEPTH of z-curve = NUMBER OF MAX SPLITS = bits used to represent single dim.  8 seems enough to represent the shape of query.
// naming is a little confusing. But, the real meaning is we have a double in unit area [0,1]. If we want to represent the space 2^16 = 65535 cells (codes),
// each dim will have 2^8=256 slices. So, we have to *scale* double value to this range [0-2^8]
// Boylece unit area yı 65535 hucre ile temsil ettik.
// Precision daha iyi olmasini istersek DEPTH= MAX_SPLIT= 32 yapabilriz. Yani unit area yı 2^64 tane hucre ile temsile edeceğiz...

//static int MAXNUMCELLINDIM = pow(2, WIDTH) - 1;
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

		double x, y;  // a point in unit area. Set the resolution :: Koordinat sistemimiz kaça kaçlık grid olacak?
		WIDTH=8; // Her eksende 8 bit kullanalım. toplam 256*256 = 65535 kod (grid cell) kullanacagiz.
		cout << "WIDTH is " << WIDTH <<endl;
		x = 0.5;
		y = 0.5;
		int xs = x * (pow(2, WIDTH) - 1); // scale to [0,255]
		int ys = y * (pow(2, WIDTH) - 1); // scale to [0,255]

		cout << "		max morcod:"<< morton2D_32_encode(255, 255)<< endl;  // 65535
		cout << "		max morcod:"<< morton2D_64_encode(255, 255)<< endl;  // 65535. evet aynı deger. 65535 tane kod kullanmak istedik. Her iki calsma uzayı da bize bunu verdi.

		uint32_t morcod = morton2D_32_encode(xs, ys); // 65536 / 4 = 16384.  code is 16383
		cout << "		morcod of center:"<< morcod<< endl;

		cout << "		morcod:"<< morton2D_32_encode(256, 256)<< endl;  // 65535*3 =196608
		cout << "		morcod:"<< morton2D_64_encode(256, 256)<< endl;  // 65535*3 =196608. evet aynı deger. 65535 tane kod kullanmak istedik. Her iki calsma uzayı da bize bunu verdi.

		// Yani uretilen morton kodlari sol alt kose pix'in morton kodu. Cell ne kadar buyuk olursa olsun aynı cunku.!!!
		// MEsela aynı uzam alanının 2^32 cell ile temsil etmek yerine 2^64 cell ile temsil edebilriiz. Yani resolution daha yüksek.

		WIDTH = 1;// Toplam 4 kod kullanacagiz. Yani bir uzam bolgesi (i.e unit area) 4 kod ile temsil ediyoruz. Resolution dusuk.
		cout << "WIDTH is " << WIDTH <<endl;
		// round kullanmalısın ki, bir sonraki slice'a bir parça geçerse onu yakalayabilelim.!
		xs = round ( 0.7 * (pow(2, WIDTH) - 1)); // scale to [0,1]
		ys = round ( 0.5 * (pow(2, WIDTH) - 1)); // scale to [0,1]
		cout << "		xs:" << xs << ", ys:" << ys << endl;
		morcod = morton2D_32_encode(xs, ys);
		cout << "		morcod of a point in unit area:"<< morcod<< endl;
		morcod = morton2D_64_encode(xs, ys);  // Gene aynı.
		cout << "		morcod of a point in unit area:"<< morcod<< endl;

		WIDTH = 2;// Toplam 16 kod kullanacagiz. Yani bir uzam bolgesi (i.e unit area) 16 kod ile temsil ediyoruz. Resolution biraz daha iyi
		cout << "WIDTH is " << WIDTH << endl;
		// round kullanmalısın ki, bir sonraki slice'a bir parça geçerse onu yakalayabilelim.!
		xs = round(0.7 * (pow(2, WIDTH) - 1)); // scale to [0,3]
		ys = round(0.5 * (pow(2, WIDTH) - 1)); // scale to [0,3]
		cout << "		xs:" << xs << ", ys:" << ys << endl;
		morcod = morton2D_32_encode(xs, ys);
		cout << "		morcod of a point in unit area:" << morcod << endl;
		morcod = morton2D_64_encode(xs, ys);  // Gene aynı.
		cout << "		morcod of a point in unit area:" << morcod << endl;

		WIDTH = 3;// Toplam 64 kod kullanacagiz. Yani bir uzam bolgesi (i.e unit area) 4 kod ile temsil ediyoruz. Resolution dusuk.
		cout << "WIDTH is " << WIDTH << endl;
		// IMPORTANT !!!!: round kullanmalısın ki, bir sonraki slice'a bir parça geçerse onu yakalayabilelim.!
		xs = round(0.7 * (pow(2, WIDTH) - 1)); // scale to [0,7]
		ys = round(0.5 * (pow(2, WIDTH) - 1)); // scale to [0,7]
		cout << "		xs:" << xs << ", ys:" << ys << endl;
		morcod = morton2D_32_encode(xs, ys);
		cout << "		morcod of a point in unit area:" << morcod << endl;
		morcod = morton2D_64_encode(xs, ys);  // Gene aynı.
		cout << "		morcod of a point in unit area:" << morcod << endl;

		// now test decoding from the last morcode = 49
		uint_fast16_t xs1,ys1;
		morton2D_32_decode(morcod,xs1,ys1);
		cout << "decoding:" << xs1 << '\t' << ys1 << endl;  // last morcod=49 ==> xs1, ys1 = 5, 4
		cout << "Cell left-bottom corner:" << (double)  xs1/pow(2, WIDTH) << '\t' <<(double) ys1/pow(2, WIDTH) << endl;  // transform it into unit area. 7/16. (7. slice coord. value)
		cout << "Cell center:" << (double)  xs1/pow(2, WIDTH) + 1/(2*pow(2, WIDTH)) << '\t' <<(double) ys1/pow(2, WIDTH) + 1/(2*pow(2, WIDTH))<< endl;  // transform it into unit area. 7/16. (7. slice coord. value)



		cout << "Generate random numbers::" <<endl;
		x = rnd.nextUniformDouble();
		y = rnd.nextUniformDouble();
		cout << x << ", " << y << endl;

		cout << morton2D_32_encode(x * (pow(2, WIDTH) - 1), y * (pow(2, WIDTH) - 1)) << endl;

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
//		return 0;
		cout <<" END OF TESTING!!!!!!!" << endl <<endl;
	}

	double *plow, *phigh;
	double xl,xh;// = xl+dx;
	double yl,yh;// = yl+dy;

	// FINDING "morton-ranges" of a range query by using QuadTree approach.
	// Use small width for understanding with debugging!!
		WIDTH=16;

	// TEST QUERY #1: (0,0) (1,1) Test the whole area.
//	xl = 0.0;
//	yl = 0.0;
//	xh = 1.0;
//	yh = 1.0;
//	plow= new double[2]{xl,yl};
//	phigh= new double[2]{xh,yh};

	// TEST QUERY #2: (0,0) (1/2^WIDTH, 1/2^WIDTH) Test the SMALLEST cell at the bottom left.
//	xl = 0*1/pow(2,WIDTH);
//	yl = 0*1/pow(2,WIDTH);
//	xh = 1*1/pow(2,WIDTH);
//	yh = 1*1/pow(2,WIDTH);
//	plow= new double[2]{xl,yl};
//	phigh= new double[2]{xh,yh};



	// TEST QUERY:  Test the random cells. Modfiy the numbers below !!  Like (0,0) (2*1/2^WIDTH, 1/2^WIDTH)
//	xl = 1-1/pow(2,WIDTH);
//	yl = 1-3*1/pow(2,WIDTH);
//	xh = 1.0;
//	yh = 1.0;
//	plow= new double[2]{xl,yl};
//	phigh= new double[2]{xh,yh};



	// RANDOM QUERY:
	do{   // generate a random query within unit area
		xl = rnd.nextUniformDouble();
		yl = rnd.nextUniformDouble();
		double dx = rnd.nextUniformDouble(0.3, 0.4);
		double dy = rnd.nextUniformDouble(0.3, 0.4);
		xh=xl+dx;
		yh=yl+dy;
	}while(xh>1.0 || yh>1.0);

	struct timeval start_time, end_time;
	double secs;
	(void)gettimeofday(&start_time, NULL);
	// TEST yuvarlama:
//	xl=0.1;
//	yl=0.1;
//	xh=0.2;
//	yh=0.3;
	xl=0.4786022785684629354818753;
	yl=0.5610612245512953677462065;
	xh=0.8401296921292444874751482;
	yh=0.9576920172544896026067818;

	cout << " Original window query::"<<endl;
	cout<< std::setprecision(25) << xl << " "<<std::setprecision(25)<< yl << " "<<std::setprecision(25)<< xh << " "<<std::setprecision(25)<< yh << " "<<endl;

	// ALIGNMENT of query corners. Approximate doubles to multiples of 1/pow(2,16). this increase the exec. time without lose of accuracy.
	// thanks to this alignment in the findMortonRanges, it never enters "else if (rq->getIntersectingArea(*subunit) > 0) { else pos==WIDTH" section..
	// xl ve yl yi sol-alt koseye; xh ve yh'yi da sag ust koseye cekmek gerek. Bu sorgunun kesiştiği en küçük hücreyi yakalamaya imkan sagliyor.
	// ((Sag ust köşe sorgu ile kesişmeyen bir hücreyi temsil ediyor. Bunun bir zararı yok. ))
	// Yuvarlama yapmasan da findMortonRanges(...) ELSE bölgesinde en küçük hücreleri zaten yakalıyor ve aynı sonuç üretiliyor.
	// Acab hangisi daha hızlı ?? Yuvarlama yaparsan 0.13 sec. Yapmazsan 0.15 sec. Yani yuvarlamada yaklaşık "10 msec" daha hızlı oluyor.
	int yuvarlama=(int)(xl/(1/pow(2,WIDTH)));
	xl=(1/pow(2,WIDTH))*yuvarlama;

	yuvarlama=(int)(yl/(1/pow(2,WIDTH)));
	yl=(1/pow(2,WIDTH))*yuvarlama;

	yuvarlama=ceil((xh/(1/pow(2,WIDTH))));  // sag-ust köşeye cekiyoruz. bu daha doğru.
//	yuvarlama=(int)(xh/(1/pow(2,WIDTH)));   // sol-alt köşeye cekiyoruz.
	xh=(1/pow(2,WIDTH))*yuvarlama;

	yuvarlama=ceil((yh/(1/pow(2,WIDTH))));
//	yuvarlama=(int)(yh/(1/pow(2,WIDTH)));
	yh=(1/pow(2,WIDTH))*yuvarlama;

	cout << " window query after yuvarlama::"<<endl;
	cout<< std::setprecision(25) << xl << " "<<std::setprecision(25)<< yl << " "<<std::setprecision(25)<< xh << " "<<std::setprecision(25)<< yh << " "<<endl;

	// Test case:  All values below are multiples of 1/pow(2,16), for example.
//	xl=0.0445404052734375;
//	xh=0.3786163330078125;
//	yl=0.100006103515625;
//	yh=0.45556640625;

	plow = new double[2]{xl,yl};
	phigh= new double[2]{xh,yh};

	cout << "Morton codes of (xl,yl) and  (xh,yh):: ";
	cout << morton2D_32_encode((uint16_t)(xl*pow(2, WIDTH)),(uint16_t)(yl*pow(2, WIDTH))) << " -- ";
	cout << morton2D_32_encode((uint16_t)(xh*pow(2, WIDTH)),(uint16_t)(yh*pow(2, WIDTH))) << endl;


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

	int pos=1;  // i.e. "split no". start using 1 bit for each dim.
	uint32_t mc=0;  // initial area is whole area, represented with mc=0
	currmin=0;
	currmax=0;

	findMortonRanges(query,unit,mc,pos);  // here traverse the quad-tree recursively
	(void)gettimeofday(&end_time, NULL);
	secs =
		(((double)end_time.tv_sec * 1000000 +
		end_time.tv_usec) -
		((double)start_time.tv_sec * 1000000 +
		start_time.tv_usec)) / 1000000;
	printf("[STAT] Finding morton codes takes %.2f seconds: ", secs);

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
//	for (auto itr = resultset.begin(); itr != resultset.end();itr++) {
//		cout << '\t' << itr->first
//				<< '\t' << itr->second << '\n';
//	}
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

			double xx1=(double)x1/pow(2, WIDTH);
			double yy1=(double)y1/pow(2, WIDTH);
			double xx2=(double)x2/pow(2, WIDTH);
			double yy2=(double)y2/pow(2, WIDTH);

//			cout << std::setprecision(10) << xx1 << '\t' << std::setprecision(10)<< yy1 << endl;
//			cout << std::setprecision(10) << xx2 << '\t' <<std::setprecision(10) << yy2 << endl;

			// convert 16-bit morton code to 32-bit codes:  2 ways: 1-transform by scaling over morton space 2- transform for each double dimension
//			uint32_t bit32mc1=( itr->first / pow(2, 2*WIDTH) ) * pow(2, 4* WIDTH);
//			uint32_t bit32mc2=( itr->second / pow(2,2*WIDTH) ) * pow(2, 4* WIDTH);
//
//			uint32_t bit32mc11= morton2D_32_encode(xx1*pow(2, 2*WIDTH), yy1*pow(2, 2*WIDTH));
//			uint32_t bit32mc22=morton2D_32_encode(xx2*pow(2, 2*WIDTH), yy2*pow(2, 2*WIDTH));
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
			// subunit'in icerisindeki (varsa) butun alt hucreleri  bulmam gerek.
			//pos = Width'e ulastıysa zaten en dipteyiz. Aksi takdirde (WIDTH-pos)*2 tane icerse alt hucreler var. (her eksende 1 split)
			int shift = (2*WIDTH - 2 * pos);
			// Su andaki quad'i shift kadar sola kaydırırsak min degerini buluyoruz.
			uint32_t mortonmin = (quad0 << shift);
			// shift ettigimiz bolgedeki bitleri 1 ile doldurursak bu quad'in en buyuk alt hucresine erismis olacagiz.
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
			if (pos < WIDTH)
				findMortonRanges(rq, subunit, quad0, pos + 1);
			else {// NEVER ENTERs HERE IFF CORNERS OF RANGE QUERY is ALIGNED TO MULTIPLES OF 1/pow(2,16). WIDTH=16
				int shift = (2*WIDTH - 2 * pos);// shift will always be ZERO here. Because we are at the end of splits, i.e. at the smallest cell.
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
			int shift = (2*WIDTH - 2 * pos);
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
			if (pos < WIDTH)
				findMortonRanges(rq, subunit, quad1, pos + 1);
			else {// NEVER ENTERs HERE IFF CORNERS OF RANGE QUERY is ALIGNED TO MULTIPLES OF 1/pow(2,16). WIDTH=16
				int shift = (2*WIDTH - 2 * pos);  // shift will always be ZERO here. Because we are at the end of splits, i.e. at the smallest cell.
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
			int shift = (2*WIDTH - 2 * pos);
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
			if (pos < WIDTH)
				findMortonRanges(rq, subunit, quad2, pos + 1);
			else { // NEVER ENTERs HERE IFF CORNERS OF RANGE QUERY is ALIGNED TO MULTIPLES OF 1/pow(2,16). WIDTH=16
				int shift = (2*WIDTH - 2 * pos);
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
			int shift = (2*WIDTH - 2 * pos);  // subunit'in icerisindeki (varsa) butun alt hucreleri  bulmam gerek. pos = Width'e ulastıysa zaten en dipteyiz. Aksi takdirde (WIDTH-pos)*2 tane icerse alt hucreler var. (her eksende 1 split)
			uint32_t mortonmin = (quad3 << shift); // Su andaki quad'i shift kadar sola kaydırırsak min degerini bulduk
			uint32_t mask = pow(2, shift) - 1;  // shift ettigimiz bolgedeki bitleri 1 ile doldurursak bu quad'in en buyu alt hucresine erismis olacagiz.
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
			if (pos < WIDTH)
				findMortonRanges(rq, subunit, quad3, pos + 1);
			else {  // NEVER ENTERs HERE IFF CORNERS OF RANGE QUERY is ALIGNED TO MULTIPLES OF 1/pow(2,16). WIDTH=16
				int shift = (2*WIDTH - 2 * pos);
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
