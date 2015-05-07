
if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi


LUMINOSITY=19.365

TEST="test"

NJETS=$1

CHANNELS="All SF OF EE MuE EMu MuMu "


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
WWTo2L2Nu_nonSkim_pow \
WWTo2L2Nu_mcnlo_nnll  \
WWTo2L2Nu_mad_nnll \
Top                \

"

#rm -rf rootfiles/${NJETS}jet

mkdir rootFiles

# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 
	
	mkdir rootFiles/${NJETS}jet/
	mkdir rootFiles/${NJETS}jet/${CHANNEL}	
	root -l -b -q "RunPROOF_test.C($LUMINOSITY,\"$TEST\",\"$SAMPLE\","$NJETS",\"$CHANNEL\")"
	mv ${SAMPLE}.root rootFiles/${NJETS}jet/${CHANNEL}
  
    done

    OUTPATH=rootFiles/${NJETS}jet/${CHANNEL}

#    Merge some MC files
    hadd -f ${OUTPATH}/WW.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu.root 
    hadd -f ${OUTPATH}/WW_pow.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu_pow.root 
    hadd -f ${OUTPATH}/WW_mc.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu_mcnlo.root
    hadd -f ${OUTPATH}/Wgamma.root ${OUTPATH}/WgammaNoStar.root ${OUTPATH}/WgammaStar.root ${OUTPATH}/GamGamWW.root 

done


