#!/usr/env bash

DATADIR="/home/jwalls/workspace/libFDR/"
WORKDIR=$DATADIR

#mincreshape -float -2 $DATADIR/data.mnc $DATADIR/data2.mnc
#$WORKDIR/forward3dfft $DATADIR/data.mnc $DATADIR/data_fft.mnc
#$WORKDIR/fftshift $DATADIR/data_fft.mnc $DATADIR/shifted_fft.mnc
#$WORKDIR/buildbwfilter $DATADIR/shifted_fft.mnc 0.9 $DATADIR/bwfilter.mnc
#$WORKDIR/buildrollofffilter $DATADIR/shifted_fft.mnc 0.3 1.0 $DATADIR/rollofffilter.mnc
#$WORKDIR/buildinversefdr $DATADIR/shifted_fft.mnc $DATADIR/normpsf_fft.mnc $DATADIR/inverse_filter.mnc
#$WORKDIR/buildpowerspectrum $DATADIR/shifted_fft.mnc $DATADIR/powerspec.mnc
#$WORKDIR/buildwienerfilter $DATADIR/powerspec.mnc $DATADIR/wienerfilter.mnc
#$WORKDIR/buildinversefdr $DATADIR/shifted_fft.mnc $DATADIR/normpsf_fft.mnc $DATADIR/inverse_filter2.mnc
#$WORKDIR/limitfdr $DATADIR/inverse_filter2.mnc 10.0 $DATADIR/limited_inverse_filter.mnc

#mincmath -mult -clobber $DATADIR/wienerfilter.mnc $DATADIR/bwfilter.mnc $DATADIR/rollofffilter.mnc $DATADIR/WBR.mnc

$WORKDIR/complexrealmultiply $DATADIR/shifted_fft.mnc $DATADIR/WBR.mnc $DATADIR/data_fft_wbr.mnc
$WORKDIR/complexmultiply $DATADIR/data_fft_wbr.mnc $DATADIR/limited_inverse_filter.mnc $DATADIR/data_fft_wbr_filt.mnc

#$WORKDIR/fftshift $DATADIR/data_fft_wbr_filt.mnc $DATADIR/data_fft_wbr_shifted.mnc

$WORKDIR/inverse3dfft $DATADIR/data_fft_wbr_filt.mnc $DATADIR/back_fdwbr.mnc

$WORKDIR/mincmagnitude $DATADIR/back_fdwbr.mnc $DATADIR/abs_bfdwbr2.mnc

