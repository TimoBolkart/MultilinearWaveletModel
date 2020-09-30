///////////////////////////////////////////////////////////////////////////////
//
//	Profile.cpp
//
//	Source file for the CProfile class
//
//	Alan Brunton 2008
//
///////////////////////////////////////////////////////////////////////////////


#include "Profile.h"


///////////////////////////////////////////////////////////////////////////////
//actions/operations
///////////////////////////////////////////////////////////////////////////////

void CProfile::computeStatistics()
{
	int i, nTrials;
	float trial;

	nTrials = m_trials.length();

	m_variance = m_mean = 0.0;
	m_min = 1e10;
	m_max = 0.0;
	if (nTrials == 0)
	{
		m_min = 0.0;
		return;
	}

	for (i = 0; i < nTrials; i++)
	{
		trial = m_trials[i];
		m_mean += trial;
		if (trial < m_min)
			m_min = trial;
		if (trial > m_max)
			m_max = trial;
	}

	m_mean /= (float) nTrials;

	for (i = 0; i < nTrials; i++)
	{
		trial = m_trials[i];
		trial -= m_mean;
		m_variance += trial * trial;
	}

	if (nTrials > 1)
		m_variance /= (float) (nTrials - 1);
}

void CProfile::reportStatistics(FILE* pfReport)
{
	_ASSERT(pfReport != NULL);

	::fprintf(pfReport, "%s:\n", m_strName.c_str());
	::fprintf(pfReport, "%i trials\n", m_trials.length());
	::fprintf(pfReport, "\ttotal= %f s\n", m_mean * (float)m_trials.length());
	::fprintf(pfReport, "\tmean = %f s\n", m_mean);
	::fprintf(pfReport, "\tvar  = %f s^2,\tstd  = %f\n", m_variance, sqrt(m_variance));
	::fprintf(pfReport, "\tmin  = %f s\n", m_min);
	::fprintf(pfReport, "\tmax  = %f s\n", m_max);
	::fprintf(pfReport, "\n");
	::fflush(pfReport);
}



