#!/bin/bash

# Forward 3D FFT
/home/johwal/packages/libFDR/forward3dfft data.mnc data_fft.mnc 
# FFTSHIFT
/home/johwal/packages/libFDR/fftshift data_fft.mnc shift_data_fft.mnc
rm data_fft.mnc
# Antialias data
/home/johwal/packages/libFDR/antialiasdata shift_data_fft.mnc aa_shift_data_fft.mnc 300
rm shift_data_fft.mnc
# Multiply by the prebuild filter
/home/johwal/packages/libFDR/complexmultiply /mnt/archive/users/johwal/fdr_BR_6x_5.mnc aa_shift_data_fft.mnc filtered_aa_shift_data_fft.mnc
rm aa_shift_data_fft.mnc
# FFTSHIFT
/home/johwal/packages/libFDR/fftshift filtered_aa_shift_data_fft.mnc shift_filtered_aa_shift_fft_data.mnc
rm filtered_aa_shift_data_fft.mnc
# Inverse 3D FFT
/home/johwal/packages/libFDR/inverse3dfft shift_filtered_aa_shift_fft_data.mnc back_shift_filtered_aa_shift_fft_data.mnc
rm shift_filtered_aa_shift_fft_data.mnc
# Magnitude
/home/johwal/packages/libFDR/mincmagnitude back_shift_filtered_aa_shift_fft_data.mnc abs_back_shift_filtered_aa_shift_fft_data.mnc 
rm back_shift_filtered_aa_shift_fft_data.mnc   