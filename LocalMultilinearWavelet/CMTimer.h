#ifndef __AB_CMTIMER_H__
#define __AB_CMTIMER_H__

class CMTimer
{
public:
CMTimer() :
	m_usePerformanceCounterMode(true),
	m_lastTimeSmall(0),
	m_currentTimeSmall(0)
{
	LARGE_INTEGER	milliSeconds;
	milliSeconds.QuadPart=(unsigned int)1000;
	m_lastTimeLarge.QuadPart=0;
	m_currentTimeLarge.QuadPart=0;
	m_elapsedMilliSeconds=0;

	if ( ! QueryPerformanceFrequency( &m_ticksPerMilliSecond ) )
	{
//		CM_LOG1(LOG_ERROR,"Does not support QueryPerformanceFrequency, returned error: %d", GetLastError() );
		m_usePerformanceCounterMode = false;
	}
	
	m_ticksPerMilliSecond.QuadPart = (__int64)m_ticksPerMilliSecond.QuadPart / (__int64)milliSeconds.QuadPart;
}

void tick()
{
	if ( m_usePerformanceCounterMode )
	{
		m_lastTimeLarge = m_currentTimeLarge;
		if ( ! QueryPerformanceCounter( &m_currentTimeLarge ) )
		{
//			CM_LOG1(LOG_ERROR,"QueryPerformanceCounter failed, reverting to timeGetTime, returned error: %d", GetLastError() );
			m_usePerformanceCounterMode = false;
			goto USE_GETTIME;
		}
		m_elapsedMilliSeconds = (int)( ( (__int64)m_currentTimeLarge.QuadPart - (__int64)m_lastTimeLarge.QuadPart ) / (__int64)m_ticksPerMilliSecond.QuadPart );
	}
	else
	{
USE_GETTIME:
		m_lastTimeSmall = m_currentTimeSmall;
		m_currentTimeSmall = timeGetTime();
		m_elapsedMilliSeconds = m_currentTimeSmall - m_lastTimeSmall;
	}
}

inline int CMTimer::getDeltaMilliSeconds()
{
	return m_elapsedMilliSeconds;
}

inline int CMTimer::getTickTime()
{
	if (m_usePerformanceCounterMode)
		return (int) ((__int64)m_currentTimeLarge.QuadPart / (__int64)m_ticksPerMilliSecond.QuadPart);
	return m_currentTimeSmall;
}

protected:
	LARGE_INTEGER	m_lastTimeLarge;
	LARGE_INTEGER	m_currentTimeLarge;
	LARGE_INTEGER	m_ticksPerMilliSecond;
	DWORD			m_lastTimeSmall;
	DWORD			m_currentTimeSmall;

	int				m_elapsedMilliSeconds;
	bool			m_usePerformanceCounterMode;
};

#endif //__AB_CMTIMER_H__




