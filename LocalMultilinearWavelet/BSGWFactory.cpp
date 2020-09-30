///////////////////////////////////////////////////////////////////////////////
//
//	BSGWFactory.cpp
//
//	Source file for the CBSGWFactory class
//
//	Alan Brunton, October 2009
//
///////////////////////////////////////////////////////////////////////////////


#include "StereoFace.h"
#include "BSGWFactory.h"


void* CBSGWFactory::createBSGW(EBSGWType eType, int nBaseWidth, int nBaseHeight, int nLevels)
{
	switch (eType)
	{
	case eBSGWFloat:
		return (void*)(new CBSplineGridWavelet<float>(nBaseWidth, nBaseHeight, nLevels, eType));
	case eBSGWDouble:
		return (void*)(new CBSplineGridWavelet<double>(nBaseWidth, nBaseHeight, nLevels, eType));
	case eBSGWFloat3:
		return (void*)(new CBSplineGridWavelet<abutil::C3Vectorf>(nBaseWidth, nBaseHeight, nLevels, eType));
	}
	
	return NULL;
}

