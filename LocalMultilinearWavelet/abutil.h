///////////////////////////////////////////////////////////////////////////////
//
//	abutil.h
//
//	Header file for the abutil group of header files
//	with additional useful includes
//
//	Alan Brunton 2009
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __ABUTIL_H__
#define __ABUTIL_H__


///////////////////////////////////////////////////////////////////////////////
//disable some warnings
///////////////////////////////////////////////////////////////////////////////

#pragma warning (disable : 4305)
#pragma warning (disable : 4244)
#pragma warning (disable : 4996)


///////////////////////////////////////////////////////////////////////////////
//includes
///////////////////////////////////////////////////////////////////////////////

//Standard library
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <cfloat>

//STL
#include <string>
#include <iostream>
#include <fstream>


typedef std::string string;


#if defined(WIN32)
#include "windows.h"
#include "minmax.h"

#define isnan(x) ((x) != (x))

#else
#endif //defined(WIN32)


//alternative definition of epsilon
#define ABUTIL_EPSILON			(1e-10)

#define __AB_CUDIFY__

namespace abutil
{

const double					g_epsilon = 1e-10;
const double					g_pi = 3.1415926535897932;

//utility code
#include "Vector.h"
#include "Quaternion.h"
#include "Matrix.h"
//#include "RandomSample.h"
//#include "Triangle.h"
#include "XArray.h"
#include "CMTimer.h"
//#include "Thread.h"
//#include "Mutex.h"
#include "Report.h"

//tokenize a string according to a set of delimiting characters
void tokenize(char* szInput, const char* cszDelim, abutil::CXArray<char*>& tokens);

};

#endif //__ABUTIL_H__


