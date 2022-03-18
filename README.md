# BotTop2PlusMinus
Convert Illumina Manifest file with Bot/Top annotation to SNP array DB compatible +/- annotation. This is only necessary if you set up a new array for the SNP array DB. 

run like:
BotTop2PlusMinus hg19.2bit manifest.csv


choose a 2 bit file build that corresponds to the manifest.csv (hg19.bit with GRCh37 or hg38.bit with GRCh38)
the manifest file comes from the genotyping lab. Load the Genomestudio project and choose export manifest. Before using the manifes run "dos2unix manifest.csv". In case you did not the following error will occur: 
"Error:Couldn't find data start"
