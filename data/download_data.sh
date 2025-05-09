#! /bin/sh

# mzTab data
wget https://raw.githubusercontent.com/HUPO-PSI/mzTab/refs/heads/master/examples/1_0-Proteomics-Release/SILAC_SQ.mzTab

# maxQuant output
wget -O proteinGroups.txt https://zenodo.org/records/3774452/files/MaxQuant_Protein_Groups.tabular?download=1

# FragPipe
wget -O data.zip --content-disposition "https://www.dropbox.com/sh/38rgczrcmdrcnm9/AADCiLYKhfouSFk-LD4gR-D0a?dl=1"
unzip data.zip combined_protein.tsv
rm data.zip

