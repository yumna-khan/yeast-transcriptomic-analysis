for sample in *.fastq.gz
do
	base=$(basename $sample .fastq.gz)
	salmon quant \
	-i combined_index \
	-l A \
	-r ${base}.fastq.gz \
    	--validateMappings \
    	-o quants/${base}_quant
done

