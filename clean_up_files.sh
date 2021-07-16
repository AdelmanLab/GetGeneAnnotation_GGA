#!/bin/bash

expname=$1
admin=$2


#make a SAF file from the BED file 
#awk '{OFS="\t";print($4,$1,$2,$3,$6)}' ${expname}_Dominant.TSS.TES.calls.bed > ${expname}_Dominant.TSS.TES.calls.saf


if [ "$admin" -eq 0 ]; then
echo "Admin mode off, removing unwanted files"

#rm TSScall files
rm ${expname}_Dominant.TSS.TES.calls_WORKING_FILE.including.unstable.txt

rm ${expname}_Dominant.TSS.TES.calls.including.unstable.bed
rm ${expname}_Dominant.TSS.TES.calls_formakeheatmap.including.unstable.txt
rm ${expname}_Dominant.TSS.calls.toLast.anno.TES.only.unstable.bed
rm ${expname}_Dominant.TSS.calls.toLast.anno.TES_formakeheatmap.only.unstable.txt


rm ${expname}_Conflict_Filtered.txt
rm ${expname}_FINAL_Dominant_obsTSS_and_TES_perGeneID_List.txt
rm ${expname}_FINAL_Dominant_obsTSS_perGeneID_List.txt
rm ${expname}_Flagged_Nondominant_obsTSS_perGeneID_List.txt
rm ${expname}_Merged.KALLISTO.ONLY.TSS.TES.calls.txt
rm ${expname}_Merged.kallisto.TSScall.TSS.TES.calls.txt
rm ${expname}_Temp_Dominant_obsTSS_perGeneID_List.txt
rm ${expname}_TSScall_detail_file
rm ${expname}_TSScall_detail_file_TSSclassify
rm ${expname}_TSScall_Output_annoTSS_List.txt
rm ${expname}_TSScall_output.bed
rm ${expname}_TSScall_Output_nuTSS_List_IGNORE.txt
rm ${expname}_TSScall_Output_obsTSS_List_Filtered_For_Dom_TSScalling.txt
rm ${expname}_TSScall_Output_obsTSS_List.txt
rm ${expname}_Weird_conflict_duplicates_obsTSS_and_TES_perGeneID_List.txt


rm -r ${expname}


#rm unnecessary Kallisto files 
rm ${expname}_TES_cluster_counts.txt 
rm ${expname}_avg_tpms_all_tx.txt 
rm ${expname}_Kallisto_called_dominant_transcript.per.gene.gtf
rm ${expname}_Kallisto_called_dominant_TSS_to_TES.bed
rm ${expname}_Kallisto_called_dominant_TSS_to_TES_clusters.txt
rm ${expname}_Kallisto_called_dominant_TSS_to_TES_cluster.transcripts.gtf
rm ${expname}_Kallisto_called_dominant_TSS_to_TES_formakeheatmap.txt
rm ${expname}_kallisto_domTSS_TES_cluster_counts.txt
rm ${expname}_avg_tpms_kallisto-filt-tx.txt
rm ${expname}_kallisto-filt.gtf


#move bash scripts (to sbatch R scripts) into logs folder
mv tsscall_cleanup_rscript.sh logs/
mv tsscall_TES_rscript.sh logs/
mv kallisto_tss_tes_rscript.sh logs/
mv gtf_filter_script.sh logs/
mv kallisto_tsscall_merge_rscript.sh logs/
mv make_tss_to_tes_transcripts_gtf_rscript.sh logs/





fi


if [ "$admin" -eq 1 ]; then
echo "Admin mode on, moving secondary files to ${expname} directory"




#rm TSScall files 
mv ${expname}_Dominant.TSS.TES.calls_WORKING_FILE.including.unstable.txt ${expname}/
mv ${expname}_Dominant.TSS.TES.calls.including.unstable.bed ${expname}/
mv ${expname}_Dominant.TSS.TES.calls_formakeheatmap.including.unstable.txt ${expname}/
mv ${expname}_Dominant.TSS.calls.toLast.anno.TES.only.unstable.bed ${expname}/
mv ${expname}_Dominant.TSS.calls.toLast.anno.TES_formakeheatmap.only.unstable.txt ${expname}/


mv ${expname}_Conflict_Filtered.txt ${expname}/
mv ${expname}_FINAL_Dominant_obsTSS_and_TES_perGeneID_List.txt ${expname}/
mv ${expname}_FINAL_Dominant_obsTSS_perGeneID_List.txt ${expname}/
mv ${expname}_Flagged_Nondominant_obsTSS_perGeneID_List.txt ${expname}/
mv ${expname}_Merged.KALLISTO.ONLY.TSS.TES.calls.txt ${expname}/
mv ${expname}_Merged.kallisto.TSScall.TSS.TES.calls.txt ${expname}/
mv ${expname}_Temp_Dominant_obsTSS_perGeneID_List.txt ${expname}/
mv ${expname}_TSScall_detail_file ${expname}/
mv ${expname}_TSScall_detail_file_TSSclassify ${expname}/
mv ${expname}_TSScall_Output_annoTSS_List.txt ${expname}/
mv ${expname}_TSScall_output.bed ${expname}/
mv ${expname}_TSScall_Output_nuTSS_List_IGNORE.txt ${expname}/
mv ${expname}_TSScall_Output_obsTSS_List_Filtered_For_Dom_TSScalling.txt ${expname}/
mv ${expname}_TSScall_Output_obsTSS_List.txt ${expname}/
mv ${expname}_Weird_conflict_duplicates_obsTSS_and_TES_perGeneID_List.txt ${expname}/

#rm unnecessary Kallisto files 
mv ${expname}_TES_cluster_counts.txt ${expname}/
mv ${expname}_avg_tpms_all_tx.txt ${expname}/
mv ${expname}_Kallisto_called_dominant_transcript.per.gene.gtf ${expname}/
mv ${expname}_Kallisto_called_dominant_TSS_to_TES.bed ${expname}/
mv ${expname}_Kallisto_called_dominant_TSS_to_TES_clusters.txt ${expname}/
mv ${expname}_Kallisto_called_dominant_TSS_to_TES_cluster.transcripts.gtf ${expname}/
mv ${expname}_Kallisto_called_dominant_TSS_to_TES_formakeheatmap.txt ${expname}/
mv ${expname}_kallisto_domTSS_TES_cluster_counts.txt ${expname}/
mv ${expname}_avg_tpms_kallisto-filt-tx.txt ${expname}/
mv ${expname}_kallisto-filt.gtf ${expname}/

#move bash scripts (to sbatch R scripts) into logs folder
mv tsscall_cleanup_rscript.sh logs/
mv tsscall_TES_rscript.sh logs/
mv kallisto_tss_tes_rscript.sh logs/
mv gtf_filter_script.sh logs/
mv kallisto_tsscall_merge_rscript.sh logs/
mv make_tss_to_tes_transcripts_gtf_rscript.sh logs/




fi





