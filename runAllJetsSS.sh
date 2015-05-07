
if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi


LUMINOSITY=19.365

TEST="test"

NJETS=$1

CHANNELS="OF"
#"All SF OF EE MuE EMu MuMu "


SAMPLES="
QCD                \ 
Top                \
WJets              \
TTJets             \
"

#rm -rf rootfiles/${NJETS}jet

mkdir rootFilesSS

# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 
	
	mkdir rootFilesSS/AllJet/
	mkdir rootFilesSS/AllJet/${CHANNEL}	
	root -l -b -q "RunPROOF_test.C($LUMINOSITY,\"$TEST\",\"$SAMPLE\","$NJETS",\"$CHANNEL\")"
	mv ${SAMPLE}.root rootFilesSS/AllJet/${CHANNEL}
  
    done

#    OUTPATH=rootFiles/AllJet/${CHANNEL}

#    Merge some MC files
#    hadd -f ${OUTPATH}/WW.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu.root 
#    hadd -f ${OUTPATH}/WW_pow.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu_pow.root 
#    hadd -f ${OUTPATH}/WW_mc.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu_mcnlo.root
#    hadd -f ${OUTPATH}/Wgamma.root ${OUTPATH}/WgammaNoStar.root ${OUTPATH}/WgammaStar.root ${OUTPATH}/GamGamWW.root 

done


