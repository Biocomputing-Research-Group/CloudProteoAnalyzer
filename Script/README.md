# Script files for identification of pFind3 and FragPipe
All scripts run on Python version 3.8.16.
## Idnetification of FragPipe
1. Copy the results of Percolator from FragPipe to working folder.
   ```
   /home/fragpipe/interact-110714_yeast_ups1_2fmol_r1.pep.xml
   ```
2. Download scipt files into working folder.
   ```
   /home/fragpipe/fragpipe2sipros.py
   /home/fragpipe/sipros_peptides_assembling.py
   /home/fragpipe/sipros_post_module.py
   ```
3. Download CloudProteoAnalyzer configuration file into working folder.
   ```
   /home/fragpipe/yeastups.cfg
   ```
4. PSMs, peptides, and proteins identification
   Use below commnad to control PSMs and peptide FDR at 1%. You can change the FDR at peptide level to make protein FDR around 1%.
   ```
   python3 pfind2sipros.py {working folder} {fdr}
   ```
   EX:
   ```
   python3 pfind2sipros.py /home/fragpipe/ 0.01
   ```
   Use below command to assemlbe peptide.
   ```
   python3 sipros_peptides_assembling.py -d {decoy label} -w {working folder} -c {cofiguration file}
   ```
   EX:
   ```
   python3 sipros_peptides_assembling.py -d DECOY_ -w /home/fragpipe/ -c /home/fragpipe/yeastups.cfg
   ```
## Idnetification of pFind3
1. Extract PSM score files from pFind3 to working folder.
   ```
   /home/pfind/yeastups.txt
   ```
2. Download scipt files into working folder.
   ```
   /home/pfind/pfind2sipros.py
   /home/pfind/sipros_peptides_assembling.py
   /home/pfind/sipros_post_module.py
   ```
3. Download CloudProteoAnalyzer configuration file into working folder.
   ```
   /home/pfind/yeastups.cfg
   ```
4. PSMs, peptides, and proteins identification
   Use below commnad to control PSMs and peptide FDR at 1%. YOu can change the FDR at peptide level to make protein FDR around 1%.
   ```
   python3 pfind2sipros.py {working folder} {fdr}
   ```
   EX:
   ```
   python3 pfind2sipros.py /home/pfind/ 0.01
   ```
   Use below command to assemlbe peptide.
   ```
   python3 sipros_peptides_assembling.py -d {decoy label} -w {working folder} -c {configuration file}
   ```
   EX:
   ```
   python3 sipros_peptides_assembling.py -d REV_ -w /home/pfind/ -c /home/pfind/yeastups.cfg
   ```
