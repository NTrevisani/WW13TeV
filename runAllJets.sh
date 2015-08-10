
if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
     echo "  "
     exit -1
fi


LUMINOSITY=0.04003

TEST="test"

NJETS=$1

CHANNELS="OF"
#"MuMu OF EE All SF OF EE MuE EMu MuMu"

PROOFMODE="kCluster"

MUONIDS="MediumIDTighterIP"
# TightID TightIDTighterIP MediumIDTighterIP"

SAMESIGN="OS" 

SAMPLES="
Data2015             \
WW50                 \
WJets50              \
WZ50                 \
ZZ50                 \
TTJets50             \
DY50                 \
VV50                 \
"
#TTbar50              \
#WW50                 \
#WJets50              \
#"
#WZ50                \
#ZZ50                \
#ZZ                  \
#WW                  \
#WJets               \

#QCD                \ 
#Top                \
#TTJets             \
#VBF                \
#WW                 \

#rm -rf rootfiles/${NJETS}jet

mkdir rootFiles

# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 
	
	for MUONID in $MUONIDS; do 
	
	    mkdir rootFiles/
	    mkdir rootFiles/${CHANNEL}	
	    mkdir rootFiles/${CHANNEL}/${MUONID}	
	    root -l -b -q "RunPROOF_test.C($LUMINOSITY,\"$TEST\",\"$SAMPLE\","$NJETS",\"$CHANNEL\",\"$PROOFMODE\",\"$SAMESIGN\",\"$MUONID\")"
	    mv ${SAMPLE}.root rootFiles/${CHANNEL}/${MUONID}
	    
	done
	
    done

done

