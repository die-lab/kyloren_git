#apply the command at the end of the analyses to remove every undesiderable lines from the quasi fasta
for i in $(ls *.txt)
do
sed 's/\#.\+//g' $i | grep -v '^$' > tmp && mv tmp $i
done
