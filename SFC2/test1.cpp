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

	cout << morton2D_64_encode(0, 0) << endl;

	Tools::Random rnd;
	return 0;
}


