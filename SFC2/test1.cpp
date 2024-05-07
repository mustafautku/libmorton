/*
 * test1.cpp
 *
 *  Created on: Aug 18, 2022
 *      Author: mustafautku
 */


#include "morton.h"  // add -I/home/mustafautku/git/libmorton/include/libmorton  "to GCC C++ compile"
#include <spatialindex/SpatialIndex.h> // add spatialindex "to GCC C++ link"

using namespace std;
using namespace libmorton;
using namespace SpatialIndex;

int main(int argc, char *argv[]) {

	cout << pow(2,16)<<endl;  // 65536
	cout << pow(2,32)<<endl;  // 4294967296

	// bu asagidakileri kullanmaliyim..
	cout << std::numeric_limits<uint16_t>::max() << endl;  // 65535
	cout<<std::numeric_limits<uint32_t>::max()<<endl;  // 4294967295

	cout << morton2D_32_encode(0, 0) << endl;
	cout << morton2D_32_encode(1, 0) << endl;
	cout << morton2D_32_encode(0, 1) << endl;
	cout << morton2D_32_encode(0, 100) << endl;
	cout << morton2D_32_encode(65535, 65535) << endl;  //4294967295
	uint_fast16_t x,y;  // for better performance this type uses 64 bit. yani unsigned long int.
	cout << sizeof(unsigned long int) << " Byte" << endl;
	morton2D_32_decode(4294967295,x,y);
	cout<< x <<","<<y<< endl;
	morton2D_32_decode(4294967296,x,y);
		cout<< x <<","<<y<< endl;

	cout << morton2D_32_encode(65535, 65536) << endl;   // seems wrong..
	cout << morton2D_32_encode(0, 1000000) << endl;  // seems wrong..

	cout << morton2D_64_encode(0, 0) << endl;
	cout << morton2D_64_encode(1, 0) << endl;
	cout << morton2D_64_encode(0, 1) << endl;
	cout << morton2D_64_encode(0, 100) << endl;




	cout<< round(0.5)<<endl;
	cout<< round(0.4) <<endl;


	Tools::Random rnd;
	return 0;
}


