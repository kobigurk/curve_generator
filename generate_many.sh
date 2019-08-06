curve=$1
residue_search_type=$2
amount=$3

for i in `seq 0 $amount`
do
  bitsize=`shuf -i 300-800 -n 1`
  sage ../without_pairing.sage $1 bitsize ${bitsize} $residue_search_type $i.curve
done
