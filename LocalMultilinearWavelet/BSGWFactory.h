///////////////////////////////////////////////////////////////////////////////
//
//	BSGWFactory.h
//
//	Header file for the CBSGWFactory class
//	responsible for creates CBSplineGridWavelet objects
//
//	Alan Brunton, October 2009
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __BSGWFACTORY_H__
#define __BSGWFACTORY_H__


#include "BSplineGridWavelet.h"

class CBSGWFactory
{
public:
	static void* createBSGW(EBSGWType eType, int nBaseWidth, int nBaseHeight, int nLevels);
};


#endif //__BSGWFACTORY_H__


