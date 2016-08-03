#!/bin/bash
# Script to to test SUPOLAR_PS picker modification



NUMBER=0000
DT=0.005
NTOUT=24000
# Seismic Unix (SU) File info
SUFILE="${NUMBER}.su"
SUFILTERED="${NUMBER}_filt.su"


# Bandpass Filter, Hz
BPLOWSTOP=0.5
BPLOWPASS=0.55
BPHIPASS=29.5
BPHISTOP=30.0
################################################################################
# Script Steps
################################################################################

# Bandpass Filter
#subfilt zerophase=0 locut=1 hicut=1 fstoplo=${BPLOWSTOP} fpasslo=${BPLOWPASS} fpasshi=${BPHIPASS} fstophi=${BPHISTOP} verbose=1 dt=0.00025 < ${SUFILE} > ${SUFILTERED}
sufilter f=${BPLOWSTOP},${BPLOWPASS},${BPHIPASS},${BPHISTOP} < ${SUFILE} > ${SUFILTERED} 
# | suwind verbose=1 key=trid min=15 max=15 | suxwigb

# Polarization Analysis
# Liberally copied from documentation for SUPOLAR

win="boxcar"                   # correlation window type
wl=3.5                           # correlation window length, sec
rlq=1.0                         # contrast factor of rectilinearity
infile=${SUFITERED}                # input 3C datafile
wpow=1                          # exponent of weighting function (RL)
dpow=1                          # exponent of directivity function
outfile="${wl}wl_"${NUMBER}     # output file (polarization filtered data)
trid=15                         # trid of component to be displayed
rl=4                            # Using Jurkevics, 1988 definition for rectilinearity
theta=3                         # Ditto for theta
kwl=100

./supolar_PS  < ${SUFILTERED} dt=${DT} wl=$wl rl=${rl} rlq=$rlq verbose=1 angle=deg theta=${theta} f1=1 rl=2 file=${outfile} win=${win}

pswigb n1=${NTOUT} < ${outfile}".pfilt" > ${outfile}"_pfilt.ps"
pswigb n1=${NTOUT} < ${outfile}".nfilt" > ${outfile}"_nfilt.ps"
pswigb n1=${NTOUT} < ${outfile}".efilt" > ${outfile}"_efilt.ps"



