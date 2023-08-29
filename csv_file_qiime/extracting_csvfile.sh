for sample in *taxa-bar-plots.qzv
do 
mv $sample ${sample%.qzv}.zip
name=${sample%_taxa-bar-plots.qzv}
if [ -d zip_dir ]
then rm -r zip_dir
mkdir zip_dir
fi
unzip -d zip_dir ${sample%.qzv}.zip
mkdir $name'_csv'
for csv_file in zip_dir/*/data/*csv
do cp $csv_file $name'_csv'/. 
cd $name'_csv'
for single_csv in *.csv
do mv $single_csv $name'_'$single_csv
done
cd ..
done
done
