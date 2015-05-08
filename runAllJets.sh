
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
VBF                \
WW                 \
QCD                \ 
Top                \
WJets              \
TTJets             \
"

#rm -rf rootfiles/${NJETS}jet

mkdir rootFiles

# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 
	
	mkdir rootFiles/AllJet/
	mkdir rootFiles/AllJet/${CHANNEL}	
	root -l -b -q "RunPROOF_test.C($LUMINOSITY,\"$TEST\",\"$SAMPLE\","$NJETS",\"$CHANNEL\")"
	mv ${SAMPLE}.root rootFiles/AllJet/${CHANNEL}
  
    done

done


