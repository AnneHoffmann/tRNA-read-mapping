# Accurate mapping of tRNA reads

This pipeline presented a best practice workflow to handle RNA-seq reads which map to multi-copy gene locations (tRNA genes, in particular) and make use of the mapped reads to call RNA-DNA differences, indicative of chemical tRNA modifications. In a preprocessing step, tRNA sequences are annotated and masked out from the genome. Additionally, pre-tRNA sequences are appended as individual chromosomes to the tRNA masked genome. To identify reads that overlap the boundaries of mature tRNAs, RNA-seq reads are mapped against the modified genome. At this step pre-tRNA reads and genomic mapped reads are filtered out. Then the remaining reads (of potential mature tRNA origin) are remapped against clustered and intronless tRNA sequences which contain CCA ends representing mature tRNAs. Only uniquely mapped reads are considered for the subsequent modification site calling. The pipeline was tested on both simulated and real life data showing a reasonable trade of between sensitivity and specificity. Detailed description and an exemplified application of the presented pipeline are reported in the manuscript [still under consideration for publication; reference will be added onces published].

The best practice workflow is available as bash script "best_practice_workflow.sh" and as Galaxy workflow "Galaxy-Workflow-Accurate_Mapping_of_tRNA_Reads.ga", repectively. Dependencies for the Galaxy workflow can be retrieved from: https://github.com/jfallmann/tRNA_ngs_mod_map_call_galaxy.


time sh best_practice_workflow.sh #TestData
real	2m10.206s, user	2m15.260s, sys	0m7.916s


For questions and help, please contact: anne.hoffmann@helmholtz-muenchen.de
