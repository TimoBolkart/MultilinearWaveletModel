///////////////////////////////////////////////////////////////////////////////
//
//	Profile.h
//
//	Header file for the CProfile class
//	used for profiling some operation
//
//	Alan Brunton 2009
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __PROFILE_H__
#define __PROFILE_H__


#include "abutil.h"
#include "XArray.h"
#include "CMTimer.h"


using namespace abutil;


class CProfile
{
protected:

	///////////////////////////////////////////////////////////////////////////
	//member variables
	///////////////////////////////////////////////////////////////////////////

	CMTimer									m_timer;
	CXArray<float>							m_trials;

	float									m_mean;
	float									m_variance;
	float									m_min;
	float									m_max;

	std::string								m_strName;


public:

	///////////////////////////////////////////////////////////////////////////
	//constructors/destructor
	///////////////////////////////////////////////////////////////////////////

	CProfile():	m_mean(0.0), m_variance(0.0), m_min(0.0), m_max(0.0) 
	{
	}


	///////////////////////////////////////////////////////////////////////////
	//accessors/modifiers
	///////////////////////////////////////////////////////////////////////////

	float getLastTrial()
	{
		if (m_trials.length() == 0)
			return 0.0;
		return m_trials.end();
	}

	void getName(std::string& str)	{ str = m_strName; }
	
	void setName(std::string str)	{ m_strName = str; }


	///////////////////////////////////////////////////////////////////////////
	//actions/operations
	///////////////////////////////////////////////////////////////////////////

	bool init()
	{
		clear();
		return true;
	}
	bool clear()
	{
		m_trials.fastInit();
		return false;
	}

	void start()				{ m_timer.tick(); }
	void stop()
	{
		m_timer.tick();
		int millis = m_timer.getDeltaMilliSeconds();
		float seconds = max(0.0, 0.001 * (float) millis);
		m_trials.pushEnd(seconds);
	}
	
	void computeStatistics();
	
	void reportStatistics(FILE* pfReport);
};


#endif //__PROFILE_H__


