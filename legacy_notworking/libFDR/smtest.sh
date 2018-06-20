#!/usr/env bash

DATADIR="/home/opt-r/Desktop/data/"
WORKDIR="/home/opt-r/packages/libFDR/"

#$WORKDIR/forward3dfft $DATADIR/data2.mnc $DATADIR/data_fft.mnc
#$WORKDIR/buildwienerfilter $DATADIR/data2.mnc $DATADIR/wienerfilter.mnc
#$WORKDIR/buildbwfilter $DATADIR/data_fft.mnc 0.9 $DATADIR/bwfilter.mnc
#$WORKDIR/buildrollofffilter $DATADIR/data_fft.mnc 0.3 1.0 $DATADIR/rollofffilter.mnc
## $WORKDIR/limitfdr $DATADIR/inverse_filter.mnc 1000.0 $WORKDIR/limited_inverse_filter.mnc
#$WORKDIR/rfftshift $DATADIR/wienerfilter.mnc $WORKDIR/ws.mnc
#$WORKDIR/rfftshift $DATADIR/bwfilter.mnc $WORKDIR/bws.mnc
#$WORKDIR/rfftshift $DATADIR/rollofffilter.mnc $WORKDIR/rs.mnc
$WORKDIR/complexrealmultiply $DATADIR/data_fft.mnc $DATADIR/ws.mnc $DATADIR/dw.mnc
$WORKDIR/complexrealmultiply $DATADIR/data_fft.mnc $DATADIR/bws.mnc $DATADIR/db.mnc
$WORKDIR/complexrealmultiply $DATADIR/data_fft.mnc $DATADIR/rs.mnc $DATADIR/dr.mnc
# $WORKDIR/complexmultiply $DATADIR/awfb.mnc $DATADIR/limited_inverse_filter.mnc $DATADIR/awfb_f.mnc
# $WORKDIR/fftshift $DATADIR/awfb_f.mnc $DATADIR/awfb_f_shifted.mnc
$WORKDIR/inverse3dfft $DATADIR/dw.mnc $DATADIR/dwback.mnc
$WORKDIR/inverse3dfft $DATADIR/db.mnc $DATADIR/dbback.mnc
$WORKDIR/inverse3dfft $DATADIR/dr.mnc $DATADIR/drback.mnc
$WORKDIR/mincmagnitude $DATADIR/dwback.mnc $DATADIR/dwfinal.mnc
$WORKDIR/mincmagnitude $DATADIR/dbback.mnc $DATADIR/dbfinal.mnc
$WORKDIR/mincmagnitude $DATADIR/drback.mnc $DATADIR/drfinal.mnc

