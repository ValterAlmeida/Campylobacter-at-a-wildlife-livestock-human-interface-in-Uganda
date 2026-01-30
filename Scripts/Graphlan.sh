Graphlan

module load Miniconda3

source activate graphlan

tail -n +2 merged_abundance_table.txt | cut -f1,3- > merged_abundance_table_reformatted.txt

export2graphlan.py --skip_rows 1 -i merged_abundance_table_reformatted.txt --tree merged_abundance.tree.txt --annotation merged_abundance.annot.txt --most_abundant 100 --abundance_threshold 0.3 --least_biomarkers 10 --annotations 5,6 --external_annotations 7 --min_clade_size 1

graphlan_annotate.py --annot merged_abundance.annot.txt merged_abundance.tree.txt merged_abundance.xml

graphlan.py --dpi 300 merged_abundance.xml merged_abundance.png --external_legends

conda deactivate graphlan