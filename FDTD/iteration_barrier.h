// SPDX-License-Identifier: GPL-3.0-or-later
// Experimental fixed-team iteration barrier; not a production default.
#ifndef OPENEMS_ITERATION_BARRIER_H
#define OPENEMS_ITERATION_BARRIER_H

#if defined(__linux__)
#include <pthread.h>
#include <stdexcept>

class IterationBarrier
{
public:
	explicit IterationBarrier(unsigned int count)
	{
		if (pthread_barrier_init(&m_barrier, NULL, count)!=0)
			throw std::runtime_error("Cannot initialize iteration barrier");
	}
	~IterationBarrier() {pthread_barrier_destroy(&m_barrier);}
	void wait()
	{
		const int result = pthread_barrier_wait(&m_barrier);
		if (result!=0 && result!=PTHREAD_BARRIER_SERIAL_THREAD)
			throw std::runtime_error("Iteration barrier failed");
	}
private:
	IterationBarrier(const IterationBarrier&) = delete;
	IterationBarrier& operator=(const IterationBarrier&) = delete;
	pthread_barrier_t m_barrier;
};
#else
#include <boost/thread/barrier.hpp>
typedef boost::barrier IterationBarrier;
#endif
#endif
