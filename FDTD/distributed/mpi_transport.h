// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_DISTRIBUTED_MPI_TRANSPORT_H
#define OPENEMS_DISTRIBUTED_MPI_TRANSPORT_H

#include "runtime.h"
#include <mpi.h>

namespace openems { namespace distributed {

// Optional adapter. MPI initialization/finalization belongs to the application.
// All methods and token progress run on the initializing MPI thread.
class MPIHaloTransport : public HaloTransport
{
public:
	MPIHaloTransport(MPI_Comm comm, const PartitionPlan&, const std::vector<uint32_t>& ranks);
	~MPIHaloTransport();
	CompletionPtr Receive(const MessageId&, HostBuffer) override;
	CompletionPtr Send(const MessageId&, HostBuffer) override;
	void Progress() override;
	void Fail(const std::string&) override;
private:
	struct Impl;
	std::shared_ptr<Impl> m_Impl;
};

class MPICollectiveContext : public CollectiveContext
{
public:
	explicit MPICollectiveContext(MPI_Comm comm);
	~MPICollectiveContext();
	double Minimum(double) override;
	double Sum(double) override;
	bool All(bool) override;
	void Fail(const std::string&) override;
private:
	MPI_Comm m_Comm;
	double Reduce(double, MPI_Op);
};

} }
#endif
