// SPDX-License-Identifier: GPL-3.0-or-later
// Test-only logical devices on ONE physical Metal GPU. Never install this plugin.
#include "FDTD/metal/gpu_backend_metal.h"
#include "gpu_build_id.h"
#include <cstdio>

namespace
{
int Count()
{
	return GPU_Backend_Metal::Plugin()->DeviceCount() ? 3 : 0;
}

bool Info(int n, GPU_DeviceInfo* info)
{
	if (n < 0 || n >= Count() || !GPU_Backend_Metal::Plugin()->DeviceInfo(0, info))
	{
		return false;
	}
	std::snprintf(
	    info->name, sizeof(info->name), "TEST logical domain %d on ONE physical Metal GPU", n);
	return true;
}

GPU_Backend* NewBackend(int n)
{
	return n >= 0 && n < Count() ? GPU_Backend_Metal::New() : nullptr;
}
}

GPU_PLUGIN_EXPORT const GPU_Plugin* openEMS_gpu_plugin()
{
	static const GPU_Plugin plugin = {
		OPENEMS_GPU_BUILD_ID, "metal-domain-test", Count, Info, NewBackend
	};
	return &plugin;
}
