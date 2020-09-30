///////////////////////////////////////////////////////////////////////////////
//
//	Report.h
//
//	Header file for the CReport class
//	used for report errors and other information both to command prompt and to file
//
//	Alan Brunton, February 2011
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __REPORT_H__
#define __REPORT_H__


#include <iostream>
#include <fstream>

#include "abutil.h"
#include "XArray.h"


using namespace abutil;


class CReport
{
protected:

	///////////////////////////////////////////////////////////////////////////
	//member variables
	///////////////////////////////////////////////////////////////////////////

	std::ofstream							m_file;


public:

	///////////////////////////////////////////////////////////////////////////
	//constructors/destructor
	///////////////////////////////////////////////////////////////////////////

	CReport(const char* szFile)
	{
		init(szFile);
	}

	~CReport()
	{
		clear();
	}


	///////////////////////////////////////////////////////////////////////////
	//accessors/modifiers
	///////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////
	//actions/operations
	///////////////////////////////////////////////////////////////////////////

	void init(const char* szFile)
	{
		_ASSERT(szFile != NULL);

		clear();

		m_file.open(szFile);
		if (!m_file.good())
			std::cout << "Error opening report file " << szFile << std::endl;
	}

	void clear()
	{
		if (m_file.is_open())
			m_file.close();
	}

	void loud(const char* sz)
	{
		m_file << sz;
		std::cout << sz;
	}
	void loud(char c)
	{
		m_file << c;
		std::cout << c;
	}
	void loud(int n)
	{
		m_file << n;
		std::cout << n;
	}
	void loud(unsigned int n)
	{
		m_file << n;
		std::cout << n;
	}
	void loud(float x)
	{
		m_file.precision(10);
		m_file << x;
		std::cout << x;
	}
	void loud(double x)
	{
		m_file.precision(15);
		m_file << x;
		std::cout << x;
	}
	void loud(const void* p)
	{
		m_file << p;
		std::cout << p;
	}

	void quiet(const char* sz)				{ m_file << sz; }
	void quiet(char c)						{ m_file << c; }
	void quiet(int n)						{ m_file << n; }
	void quiet(unsigned int n)				{ m_file << n; }
	void quiet(float x)						{ m_file.precision(10); m_file << x; }
	void quiet(double x)					{ m_file.precision(15); m_file << x; }
	void quiet(const void* p)				{ m_file << p; }

	void newline()							{ m_file << std::endl; }
	void newlineLoud()
	{
		m_file << std::endl;
		std::cout << std::endl;
	}
};


#endif //__REPORT_H__


