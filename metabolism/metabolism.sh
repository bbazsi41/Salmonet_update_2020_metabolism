#!/bin/sh

cd Salmonet2_strain_models

for i in *.json; do python3 metab_yara.py -i "$i" -o "$i"; done

wait

for j in *json.csv; do python3 merge_metabolism.py -i "$j" -o "$j" ; done

wait

mv *.json.csv_out ../Salmonet2_metab_networks

wait

rm *.json.csv

wait

cd ../Salmonet2_metab_networks
for k in *.json.csv_out ; do mv "$k" "${k%.json.csv_out}" ; done

wait

for l in * ; do python3 scripts/salmonet2_psi_mitab_converter.py -i "$l" -o results/"$l" -t "$l"  -nc 0,1 -mc 0,1  ; done

wait


#for m in *; do python3 ../scripts/salmonet2_ortholog_mapper.py -i "$m" -m ../scripts/Map-SeqNum-ID.txt -of ../scripts/OrthologousGroups.txt  --tax-id "$m" -to "$m" -o "$m"_OMA.mitab ; done
cd scripts/

rm ../results/.DS_Store

python3 locustag_to_oma.py

wait
rm ../results/.DS_Store
cut -f1-2 ../results/oma_* | sort -u >> all_met_interactions.tsv
rm ../results/.DS_Store
python locustag_to_oma2.py