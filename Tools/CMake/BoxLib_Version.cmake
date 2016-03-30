# Write the CCSEConfig.cmake file, can set as many as the following
# VERSION       = full requested version string
# VERSION_MAJOR = major version if requested, else 0
# VERSION_MINOR = minor version if requested, else 0
# VERSION_PATCH = patch version if requested, else 0
# VERSION_TWEAK = tweak version if requested, else 0
# VERSION_COUNT = number of version components, 0 to 4
SET(CCSE_VERSION_MAJOR 1)
SET(CCSE_VERSION_MINOR 3)
SET(CCSE_VERSION_PATCH 5)
SET(CCSE_VERSION_COUNT 3)
SET(CCSE_VERSION       ${CCSE_VERSION_MAJOR}.${CCSE_VERSION_MINOR}.${CCSE_VERSION_PATCH})
