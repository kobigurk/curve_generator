for i in `seq 0 10`
do
  bitsize=`shuf -i 300-800 -n 1`
  sage bls12.sage bitsize ${bitsize} $i.curve
done
