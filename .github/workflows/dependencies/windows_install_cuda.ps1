#!/usr/bin/env pwsh
#
# Copyright 2024 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

# fail on first error
$ErrorActionPreference = "Stop"

# skip if already installed
if (Test-Path "C:\Program Files\NVIDIA GPU Computing Toolkit\") { exit }

# configure
$cuda_version_long = '12.4.0'
$cuda_version_short = '12.4'

# download
New-item -ItemType directory -Name cuda
Set-Location -Path cuda -PassThru
Invoke-WebRequest -Uri "https://developer.download.nvidia.com/compute/cuda/${cuda_version_long}/network_installers/cuda_${cuda_version_long}_windows_network.exe" -OutFile "cuda_install.exe"

# install
$CudaFeatures = " nvcc_${cuda_version_short} " + `
" cublas_${cuda_version_short} " + `
" cublas_dev_${cuda_version_short} " + `
" cuda_profiler_api_${cuda_version_short} " + `
" cudart_${cuda_version_short} " + `
" cufft_${cuda_version_short} " + `
" cufft_dev_${cuda_version_short} " + `
" cuobjdump_${cuda_version_short} " + `
" cupti_${cuda_version_short} " + `
" curand_${cuda_version_short} " + `
" curand_dev_${cuda_version_short} " + `
" cusolver_${cuda_version_short} " + `
" cusolver_dev_${cuda_version_short} " + `
" cusparse_${cuda_version_short} " + `
" cusparse_dev_${cuda_version_short} " + `
" cuxxfilt_${cuda_version_short} " + `
" npp_${cuda_version_short} " + `
" npp_dev_${cuda_version_short} " + `
" nvdisasm_${cuda_version_short} " + `
" nvjitlink_${cuda_version_short} " + `
" nvjpeg_${cuda_version_short} " + `
" nvjpeg_dev_${cuda_version_short} " + `
" nvml_dev_${cuda_version_short} " + `
" nvprof_${cuda_version_short} " + `
" nvprune_${cuda_version_short} " + `
" nvrtc_${cuda_version_short} " + `
" nvrtc_dev_${cuda_version_short} " + `
" nvtx_${cuda_version_short} " + `
" sanitizer_${cuda_version_short} " + `
" thrust_${cuda_version_short} " + `
" visual_studio_integration_${cuda_version_short} "

$SilentFlag = '-s '
Start-Process -FilePath '.\cuda_install.exe' -ArgumentList @($SilentFlag + $CudaFeatures) -Wait -NoNewWindow

# cleanup
Remove-Item cuda_install.exe

# GitHub actions env updates
if (Test-Path "env:GITHUB_ACTIONS") {
    echo "Adding CUDA to CUDA_PATH and PATH"
    $env:CUDA_PATH = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v${cuda_version_short}"
    echo "CUDA_PATH=$env:CUDA_PATH" | Out-File -FilePath $env:GITHUB_ENV -Encoding utf8 -Append
    echo "$env:CUDA_PATH/bin" | Out-File -FilePath $env:GITHUB_PATH -Encoding utf8 -Append
}
