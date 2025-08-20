# EBV

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ebv_strands_match.csv \
  --out ../results/EBV/ebv_act_R.csv \
  --treat-hits-out ../results/EBV/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/EBV/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/EBV/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/EBV/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=NC_007605.1 \
     --outdir=../results/EBV/ \
     --prefix=EBV

Rscript 2_plot_heatmap.r \
  --matrix ../results/EBV/ebv_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/EBV/EBV_nearby_genes.bed \
  --genes-anno ../results/EBV/gene_annotations.tsv \
  --out ../results/EBV/EBV_mpra_heatmap.png \
  --title "MPRA Differential Activity – EBV" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/EBV/ebv_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/EBV/EBV_nearby_genes.bed \
  --genes-anno ../results/EBV/gene_annotations.tsv \
  --out ../results/EBV/EBV_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – EBV" \
  --kb-step=25000 --combined



# HCMV

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/hcmv_strands_match.csv \
  --out ../results/HCMV/hcmv_act_R.csv \
  --treat-hits-out ../results/HCMV/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/HCMV/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/HCMV/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/HCMV/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=FJ616285.1 \
     --outdir=../results/HCMV/ \
     --prefix=HCMV

Rscript 2_plot_heatmap.r \
  --matrix ../results/HCMV/hcmv_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/HCMV/HCMV_nearby_genes.bed \
  --genes-anno ../results/HCMV/gene_annotations.tsv \
  --out ../results/HCMV/HCMV_mpra_heatmap.png \
  --title "MPRA Differential Activity – HCMV" \
  --kb-step=25000

Rscript 2_plot_heatmap.r \
  --matrix ../results/HCMV/hcmv_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/HCMV/HCMV_nearby_genes.bed \
  --genes-anno ../results/HCMV/gene_annotations.tsv \
  --out ../results/HCMV/HCMV_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – HCMV" \
  --kb-step=25000 --combined



# KHSV

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/khsv_strands_match.csv \
  --out ../results/KHSV/khsv_act_R.csv \
  --treat-hits-out ../results/KHSV/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/KHSV/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/KHSV/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/KHSV/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=MZ712172.1 \
     --outdir=../results/KHSV/ \
     --prefix=KHSV

Rscript 2_plot_heatmap.r \
  --matrix ../results/KHSV/khsv_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/KHSV/KHSV_nearby_genes.bed \
  --genes-anno ../results/KHSV/gene_annotations.tsv \
  --out ../results/KHSV/KHSV_mpra_heatmap.png \
  --title "MPRA Differential Activity – KHSV" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/KHSV/khsv_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/KHSV/KHSV_nearby_genes.bed \
  --genes-anno ../results/KHSV/gene_annotations.tsv \
  --out ../results/KHSV/KHSV_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – KHSV" \
  --kb-step=25000 --combined



# VZV

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/vzv_strands_match.csv \
  --out ../results/VZV/vzv_act_R.csv \
  --treat-hits-out ../results/VZV/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/VZV/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/VZV/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/VZV/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=KU926311.1 \
     --outdir=../results/VZV/ \
     --prefix=VZV

Rscript 2_plot_heatmap.r \
  --matrix ../results/VZV/vzv_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/VZV/VZV_nearby_genes.bed \
  --genes-anno ../results/VZV/gene_annotations.tsv \
  --out ../results/VZV/VZV_mpra_heatmap.png \
  --title "MPRA Differential Activity – VZV" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/VZV/vzv_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/VZV/VZV_nearby_genes.bed \
  --genes-anno ../results/VZV/gene_annotations.tsv \
  --out ../results/VZV/VZV_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – VZV" \
  --kb-step=25000 --combined


# HSV1
python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/hsv1_strands_match.csv \
  --out ../results/HSV1/hsv1_act_R.csv \
  --treat-hits-out ../results/HSV1/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/HSV1/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/HSV1/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/HSV1/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=KT899744.1 \
     --outdir=../results/HSV1/ \
     --prefix=HSV1

Rscript 2_plot_heatmap.r \
  --matrix ../results/HSV1/hsv1_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/HSV1/HSV1_nearby_genes.bed \
  --genes-anno ../results/HSV1/gene_annotations.tsv \
  --out ../results/HSV1/HSV1_mpra_heatmap.png \
  --title "MPRA Differential Activity – HSV1" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/HSV1/hsv1_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/HSV1/HSV1_nearby_genes.bed \
  --genes-anno ../results/HSV1/gene_annotations.tsv \
  --out ../results/HSV1/HSV1_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – HSV1 KOS (combined)" \
  --kb-step=25000 --combined



# HIV1 Rejo

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/hiv1_rejo_strands_match.csv \
  --out ../results/HIV1_rejo/hiv1_rejo_act_R.csv \
  --treat-hits-out ../results/HIV1_rejo/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/HIV1_rejo/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/HIV1_rejo/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/HIV1_rejo/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=JN944911 \
     --outdir=../results/HIV1_rejo/ \
     --prefix=HIV1_rejo

Rscript 2_plot_heatmap.r \
  --matrix ../results/HIV1_rejo/hiv1_rejo_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/HIV1_rejo/HIV1_rejo_nearby_genes.bed \
  --out ../results/HIV1_rejo/HIV1_rejo_mpra_heatmap.png \
  --title "MPRA Differential Activity – HIV1 Rejo" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/HIV1_rejo/hiv1_rejo_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/HIV1_rejo/HIV1_rejo_nearby_genes.bed \
  --out ../results/HIV1_rejo/HIV1_rejo_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – HIV1 Rejo" \
  --kb-step=25000 --combined


# Adenoviruses

ad1, ad3, ad4, ad5, ad11, ad14, ad37

# ad4

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ad4_strands_match.csv \
  --out ../results/AD4/ad4_act_R.csv \
  --treat-hits-out ../results/AD4/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/AD4/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/AD4/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/AD4/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=KX384949.1 \
     --outdir=../results/AD4/ \
     --prefix=AD4

Rscript 2_plot_heatmap.r \
  --matrix ../results/AD4/ad4_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD4/AD4_nearby_genes.bed \
  --out ../results/AD4/AD4_mpra_heatmap.png \
  --title "MPRA Differential Activity – AD4" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/AD4/ad4_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD4/AD4_nearby_genes.bed \
  --out ../results/AD4/AD4_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – AD4 (combined)" \
  --kb-step=25000 --combined


# AD11

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ad11_strands_match.csv \
  --out ../results/AD11/ad11_act_R.csv \
  --treat-hits-out ../results/AD11/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/AD11/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/AD11/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/AD11/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=AF532578.1 \
     --outdir=../results/AD11/ \
     --prefix=AD11

Rscript 2_plot_heatmap.r \
  --matrix ../results/AD11/ad11_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD11/AD11_nearby_genes.bed \
  --out ../results/AD11/AD11_mpra_heatmap.png \
  --title "MPRA Differential Activity – AD11" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/AD11/ad11_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD11/AD11_nearby_genes.bed \
  --out ../results/AD11/AD11_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – AD11 (combined)" \
  --kb-step=25000 --combined


# AD3

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ad3_strands_match.csv \
  --out ../results/AD3/ad3_act_R.csv \
  --treat-hits-out ../results/AD3/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/AD3/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/AD3/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/AD3/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=AY599834.1 \
     --outdir=../results/AD3/ \
     --prefix=AD3

Rscript 2_plot_heatmap.r \
  --matrix ../results/AD3/ad3_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD3/AD3_nearby_genes.bed \
  --out ../results/AD3/AD3_mpra_heatmap.png \
  --title "MPRA Differential Activity – AD3" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/AD3/ad3_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD3/AD3_nearby_genes.bed \
  --out ../results/AD3/AD3_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – AD3 (combined)" \
  --kb-step=25000 --combined


# ad5


python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ad5_strands_match.csv \
  --out ../results/AD5/ad5_act_R.csv \
  --treat-hits-out ../results/AD5/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/AD5/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/AD5/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/AD5/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=AY601635.1\
     --outdir=../results/AD5/ \
     --prefix=AD5

Rscript 2_plot_heatmap.r \
  --matrix ../results/AD5/ad5_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD5/AD5_nearby_genes.bed \
  --out ../results/AD5/AD5_mpra_heatmap.png \
  --title "MPRA Differential Activity – AD5" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/AD5/ad5_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD5/AD5_nearby_genes.bed \
  --out ../results/AD5/AD5_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – AD5 (combined)" \
  --kb-step=25000 --combined

# ad7


python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ad7_strands_match.csv \
  --out ../results/AD7/ad7_act_R.csv \
  --treat-hits-out ../results/AD7/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/AD7/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/AD7/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/AD7/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=AY594255\
     --outdir=../results/AD7/ \
     --prefix=AD7

Rscript 2_plot_heatmap.r \
  --matrix ../results/AD7/ad7_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD7/AD7_nearby_genes.bed \
  --out ../results/AD7/AD7_mpra_heatmap.png \
  --title "MPRA Differential Activity – AD7" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/AD7/ad7_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD7/AD7_nearby_genes.bed \
  --out ../results/AD7/AD7_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – AD7 (combined)" \
  --kb-step=25000 --combined



# ad 1
python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ad1_strands_match.csv \
  --out ../results/AD1/ad1_act_R.csv \
  --treat-hits-out ../results/AD1/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/AD1/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/AD1/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/AD1/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=AC_000017.1\
     --outdir=../results/AD1/ \
     --prefix=AD1

Rscript 2_plot_heatmap.r \
  --matrix ../results/AD1/ad1_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD1/AD1_nearby_genes.bed \
  --out ../results/AD1/AD1_mpra_heatmap.png \
  --title "MPRA Differential Activity – AD1" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/AD1/ad1_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD1/AD1_nearby_genes.bed \
  --out ../results/AD1/AD1_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – AD1 (combined)" \
  --kb-step=25000 --combined



# AD14
python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ad14_strands_match.csv \
  --out ../results/AD14/ad14_act_R.csv \
  --treat-hits-out ../results/AD14/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/AD14/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/AD14/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/AD14/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=AY803294.1\
     --outdir=../results/AD14/ \
     --prefix=AD14

Rscript 2_plot_heatmap.r \
  --matrix ../results/AD14/ad14_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD14/AD14_nearby_genes.bed \
  --out ../results/AD14/AD14_mpra_heatmap.png \
  --title "MPRA Differential Activity – AD14" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/AD14/ad14_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD14/AD14_nearby_genes.bed \
  --out ../results/AD14/AD14_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – AD14 (combined)" \
  --kb-step=25000 --combined


# AD37

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ad37_strands_match.csv \
  --out ../results/AD37/ad37_act_R.csv \
  --treat-hits-out ../results/AD37/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/AD37/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/AD37/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/AD37/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=AB448775.1 \
     --outdir=../results/AD37/ \
     --prefix=AD37

Rscript 2_plot_heatmap.r \
  --matrix ../results/AD37/ad37_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD37/AD37_nearby_genes.bed \
  --out ../results/AD37/AD37_mpra_heatmap.png \
  --title "MPRA Differential Activity – AD37" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/AD37/ad37_act_R.csv \
  --conditions DMSO,H2O,HeatCtrl,DOX,ABT,Dex,H2O2,Heat,MET,MG,Oligo \
  --controls DMSO,H2O,HeatCtrl \
  --order DMSO,ABT,Dex,DOX,MG,Oligo,HeatCtrl,Heat,H2O,H2O2,MET \
  --genes-bed ../results/AD37/AD37_nearby_genes.bed \
  --out ../results/AD37/AD37_mpra_heatmap_combined.png \
  --title "MPRA Differential Activity – AD3 (combined)" \
  --kb-step=25000 --combined






# Full

# EBV

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/ebv_strands_match.csv \
  --out ../results/EBV_BB/ebv_act_R.csv \
  --treat-hits-out ../results/EBV_BB/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/EBV_BB/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/EBV_BB/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/EBV_BB/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=NC_007605.1 \
     --outdir=../results/EBV_BB/ \
     --prefix=EBV

Rscript 2_plot_heatmap.r \
  --matrix ../results/EBV_BB/ebv_act_R.csv \
  --conditions H2O,H2O2,Oligo,TNF,INFG \
  --controls H2O \
  --order H2O,H2O2,Oligo,TNF,INFG \
  --genes-bed ../results/EBV_BB/EBV_nearby_genes.bed \
  --genes-anno ../results/EBV_BB/gene_annotations.tsv \
  --out ../results/EBV_BB/EBV_mpra_heatmap.pdf \
  --title "MPRA Differential Activity – EBV" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/EBV_BB/ebv_act_R.csv \
  --conditions H2O,H2O2,Oligo,TNF,INFG \
  --controls H2O \
  --order H2O,H2O2,Oligo,TNF,INFG \
  --genes-bed ../results/EBV_BB/EBV_nearby_genes.bed \
  --genes-anno ../results/EBV_BB/gene_annotations.tsv \
  --out ../results/EBV_BB/EBV_mpra_heatmap_combined.pdf \
  --title "MPRA Differential Activity – EBV" \
  --kb-step=25000 --combined



# HCMV

python 0_preprocess_tiles.py \
  --controls-dir ../data/control_activity \
  --treatments-dir ../data/treatments_effect \
  --map ../data/strands/hcmv_strands_match.csv \
  --out ../results/HCMV_BB/hcmv_act_R.csv \
  --treat-hits-out ../results/HCMV_BB/treatment_hits.csv \
  --control-lfc 2 \
  --treat-lfc 0.585 \
  --dna-threshold 50 \
  --verbose

python 0_pull_seq.py \
  --ids ../results/HCMV_BB/treatment_hits.csv \
  --seqs ../data/all_tiles_seqs.tsv \
  --out ../results/HCMV_BB/significant_tiles_sequences.csv


Rscript 1_generate_bed.r \
     --excel=../results/HCMV_BB/significant_tiles_sequences.csv \
     --seq-col=sequence \
     --id-col=ID \
     --accession=FJ616285.1 \
     --outdir=../results/HCMV_BB/ \
     --prefix=HCMV

Rscript 2_plot_heatmap.r \
  --matrix ../results/HCMV_BB/hcmv_act_R.csv \
  --conditions H2O,H2O2,Oligo,TNF,INFG \
  --controls H2O \
  --order H2O,H2O2,Oligo,TNF,INFG \
  --genes-bed ../results/HCMV_BB/HCMV_nearby_genes.bed \
  --genes-anno ../results/HCMV_BB/gene_annotations.tsv \
  --out ../results/HCMV_BB/HCMV_mpra_heatmap.pdf \
  --title "MPRA Differential Activity – HCMV" \
  --kb-step=25000


Rscript 2_plot_heatmap.r \
  --matrix ../results/HCMV_BB/hcmv_act_R.csv \
  --conditions H2O,H2O2,Oligo,TNF,INFG \
  --controls H2O \
  --order H2O,H2O2,Oligo,TNF,INFG \
  --genes-bed ../results/HCMV_BB/HCMV_nearby_genes.bed \
  --genes-anno ../results/HCMV_BB/gene_annotations.tsv \
  --out ../results/HCMV_BB/HCMV_mpra_heatmap_combined.pdf \
  --title "MPRA Differential Activity – HCMV" \
  --kb-step=25000 --combined