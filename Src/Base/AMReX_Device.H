
#ifndef BL_DEVICE_H
#define BL_DEVICE_H

#ifdef CUDA
extern "C" {
    void initialize_cuda();
}
#endif

namespace amrex {

class Device
{

public:

    static void beginDeviceLaunchRegion() { in_device_launch_region = true; }

    static void endDeviceLaunchRegion() { in_device_launch_region = false; }

    static void setDeviceLaunchRegion(bool r) { in_device_launch_region = r; }

    static bool inDeviceLaunchRegion() {
#ifdef CUDA
      cudaDeviceSynchronize(); 
#endif
      return in_device_launch_region; 
    }

#ifdef CUDA
    static int cudaDeviceId() { return cuda_device_id; }
#endif

private:

    static bool in_device_launch_region;
    static int cuda_device_id;

};

}

#endif
