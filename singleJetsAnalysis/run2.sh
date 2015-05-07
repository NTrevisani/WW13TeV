if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi


LUMINOSITY=19.468

TEST="test"

NJETS=$1

CHANNELS="All SF OF EE MuE EMu MuM "


SAMPLES="          \
DataRun2012_Total  \
ggWWto2L           \
WWTo2L2Nu          \
WWTo2L2Nu_pow      \
WWTo2L2Nu_mcnlo    \
GamGamWW           \
WZ                 \
ZZ                 \
TTbar              \
TW                 \
DY                 \
WgammaStar         \
WgammaNoStar       \
DYtautau           \
Zgamma             \
HWW125             \
VVV                \
WJets              \
WJetsFakes_Total   \
WWTo2L2Nu_pow_nnll \
WW_pow_nnll        \
"

for CHANNEL in $CHANNELS; do

    OUTPATH=rootfiles/${NJETS}jet/${CHANNEL}

#    Merge some MC files
    hadd -f ${OUTPATH}/WW.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu.root 
    hadd -f ${OUTPATH}/WW_pow.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu_pow.root 
    hadd -f ${OUTPATH}/WW_mc.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu_mcnlo.root
    hadd -f ${OUTPATH}/Wgamma.root ${OUTPATH}/WgammaNoStar.root ${OUTPATH}/WgammaStar.root ${OUTPATH}/GamGamWW.root 

done
