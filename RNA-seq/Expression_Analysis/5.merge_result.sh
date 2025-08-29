ls ./5.Quantification/*.count > genes.quant_files.txt

perl /home/sjb/workspace/rna-seq/script/abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files genes.quant_files.txt --out_prefix ./6.Merge_result/genes
