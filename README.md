# FunFun: ITS based gene content predictor tool
<img src="https://user-images.githubusercontent.com/53526550/179017255-7efbfcb2-cee6-427a-9d43-0e2a570b2d45.png" width="200" height="280" align="left">
FunFun (Fungal Functional predictor) is the tool allows us to evaluate the functional content of an individual fungus or mycobiome based on ITS-amplicon sequencing data. 

## **Dependencies:**
- numpy
- pandas
- biopython=1.77
- plotly
- scipy
## **Installation:**
FunFun can be installed via pip:
```
conda create -n funfun_env -c bioconda
conda activate funfun_env
pip install funfun
```
Or you can clone git repository:
```
git clone https://github.com/DanilKrivonos/FunFun
python funfun --help
```
## **Ussages example:**

FunFun use three different types of ITS region input: ITS1, ITS2 and ITS1+5.8S+ITS2 (Full size ITS cluster). Also you can change model parameters: K neighbours and ε neighborhood limit (decreasing the parameters will increase the prediction accuracy). The output file is tsv tabl, when first column is KEGG orthology group and subsequent correspond to predicted gene content vectors.

```
funfun -ITS example/its2.fasta -type its2
```
In the result of FunFun you take comfortable bar plot:

<img src="https://user-images.githubusercontent.com/53526550/177324645-b297fc42-a774-4883-9cbf-8810cccea0a7.png" width="1000" height="600" align="left">

## Parameters
```
usage: FunFun.py [-h] [-ITS ITS] [-type TYPE] [-out OUT] [-K K] [-e E]

optional arguments:
  -h, --help  show this help message and exit

Base arguments:
  -ITS ITS    Sequences of ITS in fasta format.
  -type TYPE  Type of its region: its1, its2, concatenate.
  -out OUT   Output directory (default ./FunFun_output)

Algorithm Parameter:
  -K K        K nearest neighbors. Default k=10.
  -e E        Epsilon distance for nearest neighbor. 0 <= ε <= 1. Default ε=0.5.
```
## Reference
Will be add later ...
 
