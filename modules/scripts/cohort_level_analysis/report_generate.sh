for i in {1..9}; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-$i/align/mapping.csv analysis/align/mapping.csv_$i; done

#for i in {1..9}; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-$i/align/run_info.txt analysis/align/run_info.txt_$i; done

for i in {1..9}; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-$i/ceas/dhs.csv analysis/ceas/dhs.csv_$i; done

for i in {1..9}; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-$i/ceas/meta.csv analysis/ceas/meta.csv_$i; done

for i in {1..9}; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-$i/contam/contamination.csv analysis/contam/contamination.csv_$i; done

for i in {1..9}; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-$i/frips/frips.csv analysis/frips/frips.csv_$i; done

for i in {1..9}; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-$i/frips/pbc.csv analysis/frips/pbc.csv_$i; done

for i in {1..9}; do echo gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-$i/peaks/peakStats.csv analysis/peaks/peakStats.csv_$i; done


gsutil ls gs://09152021_stanford_atac_gd2car_tommy/analysis-*/conserv/CA44*/CA44*_conserv_thumb.png | cut -d / -f 6 > samples_rep.txt

cat samples.txt | while read i ; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-*/frag/$i/$i\_\frags.txt analysis/frag/$i/$i\_\frags.txt ; done

cat samples_rep.txt | while read i ; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-*/conserv/$i/$i\_\conserv_thumb.png analysis/conserv/$i/$i\_\conserv_thumb.png ; done

cat samples_rep.txt | while read i ; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-*/peaks/$i/$i\_\treat_pileup.bw analysis/peaks/$i/$i\_\treat_pileup.bw ; done

cat samples_rep.txt | while read i ; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-*/peaks/$i/$i\_\sorted_peaks.bed analysis/peaks/$i/$i\_\sorted_peaks.bed ; done

cat samples_rep.txt | while read i ; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-*/peaks/$i/$i\_\sorted_summits.bed analysis/peaks/$i/$i\_\sorted_summits.bed ; done

cat samples_rep.txt | while read i ; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-*/peaks/$i/$i\_\sorted_peaks.narrowPeak analysis/peaks/$i/$i\_\sorted_peaks.narrowPeak ; done

cat samples.txt | while read i ; do gsutil -m cp gs://09152021_stanford_atac_gd2car_tommy/analysis-*/align/$i/$i\.sorted.bam analysis/align/$i/$i\.sorted.bam ; done


head -1 mapping.csv_1 > mapping.csv
tail -n +2 -q mapping.csv_* >> mapping.csv


head -1 contamination.csv_1 > contamination.csv
tail -n +2 -q contamination.csv_* >> contamination.csv
