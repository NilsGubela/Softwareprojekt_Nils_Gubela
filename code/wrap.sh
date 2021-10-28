SEQ="CUGCGGCUUUGGCUCUAGCC"
echo $SEQ | DrTransformer --logfile

output="Next"
while [ "$output" = "Next" ]
do
	tag=$( tail -n 1 base.txt )
	echo $SEQ | DrTransformer --logfile $tag	
	output=$(python3 adaptive_step.py)
done
