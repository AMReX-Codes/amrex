
#ifndef BL_DEVICE_H
#define BL_DEVICE_H

namespace amrex {

class Device
{

public:

    static void beginDeviceLaunchRegion() { in_device_launch_region = true; }

    static void endDeviceLaunchRegion() { in_device_launch_region = false; }

    static void setDeviceLaunchRegion(bool r) { in_device_launch_region = r; }

    static bool inDeviceLaunchRegion() { return in_device_launch_region; }

    static bool in_device_launch_region;

};

}

#endif
