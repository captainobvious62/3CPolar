#!/bin/bash
# Script to to test SUPOLAR_PS picker modification



NUMBER=121855
DT=0.00025
NTOUT=64080
# Seismic Unix (SU) File info
SUFILE="${NUMBER}.su"
SUFILTERED="${NUMBER}_filt.su"


# Bandpass Filter, Hz
BPLOWSTOP=0.1
BPLOWPASS=0.5
BPHIPASS=800
BPHISTOP=820
################################################################################
# Script Steps
################################################################################

# Bandpass Filter
subfilt zerophase=0 locut=1 hicut=1 fstoplo=${BPLOWSTOP} fpasslo=${BPLOWPASS} fpasshi=${BPHIPASS} fstophi=${BPHISTOP} verbose=1 dt=0.00025 < ${SUFILE} > ${SUFILTERED}
# | suwind verbose=1 key=trid min=15 max=15 | suxwigb

# Polarization Analysis
# Liberally copied from documentation for SUPOLAR

win="hanning"                   # correlation window type
wl=0.75                         # correlation window length, sec
rlq=1.0                         # contrast factor of rectilinearity
infile=${SUFILE}                # input 3C datafile
wpow=1                          # exponent of weighting function (RL)
dpow=1                          # exponent of directivity function
outfile="${wl}wl_"${NUMBER}     # output file (polarization filtered data)
trid=15                         # trid of component to be displayed
rl=2                            # Using Jurkevics, 1988 definition for rectilinearity
theta=3                         # Ditto for theta
kwl=100

./supolar_PS  < ${infile} dt=${DT} wl=$wl rl=${rl} rlq=$rlq verbose=1 angle=deg theta=${theta} f1=1 rl=2 file=${outfile} win=${win}

pswigb n1=${NTOUT} < ${outfile}".pfilt" > ${outfile}"_pfilt.ps"
pswigb n1=${NTOUT} < ${outfile}".nfilt" > ${outfile}"_nfilt.ps"
pswigb n1=${NTOUT} < ${outfile}".efilt" > ${outfile}"_efilt.ps"



