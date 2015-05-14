
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
#VBF                \
#WW                 \

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

done


