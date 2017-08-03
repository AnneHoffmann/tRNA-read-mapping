# Accurate mapping of tRNA reads

This pipeline presented a best practice workflow to handle RNA-seq reads which map to multi-copy gene locations (tRNA genes, in particular) and make use of the mapped reads to call RNA-DNA differences, indicative of chemical tRNA modifications. In a preprocessing step, tRNA sequences are annotated and masked out from the genome. Additionally, pre-tRNA sequences are appended as individual chromosomes to the tRNA masked genome. To identify reads that overlap the boundaries of mature tRNAs, RNA-seq reads are mapped against the modified genome. At this step pre-tRNA reads and genomic mapped reads are filtered out. Then the remaining reads (of potential mature tRNA origin) are remapped against clustered and intronless tRNA sequences which contain CCA ends representing mature tRNAs. Only uniquely mapped reads are considered for the subsequent modification site calling. The pipeline was tested on both simulated and real life data showing a reasonable trade of between sensitivity and specificity. Detailed description and an exemplified application of the presented pipeline are reported in the manuscript [still under consideration for publication; reference will be added onces published].


time bash best_practice_workflow.sh 

real	1m24.750s
user	2m29.712s
sys	0m10.509s
