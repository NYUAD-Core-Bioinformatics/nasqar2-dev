#wget wget --no-check-certificate 'https://drive.google.com/u/1/uc?id=1dUE8CquO56CX8RgDnrCpU1D24nLV3UYJ&export=download' -O data.zip


#export fileid=1tpNYzDrf1KnIkElEMWNTdlb2JYVN0HfJ
export filename=data.zip

## WGET ##
#wget --save-cookies cookies.txt 'https://docs.google.com/uc?export=download&id='$fileid -O- \
#     | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p' > confirm.txt

#wget --load-cookies cookies.txt -O $filename \
#     'https://docs.google.com/uc?export=download&id='$fileid'&confirm='$(<confirm.txt)

#rm cookies.txt confirm.txt > /dev/null

unzip data.zip

cp -r data/dada2Shiny/www  dada2Shiny/src/
cp -r data/deseq2shiny/www  deseq2shiny/src/
cp -r data/ClusterProfShinyORA/www  ClusterProfShinyORA/src/
cp -r data/ClusterProfShinyGSEA/www  ClusterProfShinyGSEA/src/
cp -r data/DEApp/www  DEApp/src/
cp -r data/GeneCountMerger/www GeneCountMerger/src/
cp -r  data/NASQAR/tsar_nasqar/www      tsar_nasqar/src/
cp -r  data/NASQAR/mergeFPKMs/www      mergeFPKMs/src/
cp -r  data/ATACseqQCShniy/www      ATACseqQCShniy/src/