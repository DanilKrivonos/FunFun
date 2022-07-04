FunFun (Fungal Functional predictor)is the instrument allows us to evaluate the functional content of an individual fungus or mycobiome based on ITS-amplicon sequencing data. 
<img src="https://user-images.githubusercontent.com/53526550/177148554-e6db1eab-e4e2-4e44-90a0-66f3ba8fbb47.png" width="200" height="280" align="left">
## **Dependencies:**
- numpy
- tqdm
- pandas
- biopython=1.77
- plotly
- scipy
## **Installation:**

```
conda create -n funfun_env -c bioconda
conda activate funfun_env
pip install funfun
```

## **Ussages example:**

FunFun use three different types of ITS region input: ITS1, ITS2 and 18S+ITS1+5.8S+ITS2+28S (Full size ITS cluster). Also you can change model parameters: K neighbours and ε neighborhood limit. The output file is tsv tabl, when first column is orthology group and second is representation.

```
funfun -ITS example/its2.fasta -type its2
```
### Parameters
```
usage: FunFun.py [-h] [-ITS ITS] [-type TYPE] [-out OUT] [-K K] [-e E]

optional arguments:
  -h, --help  show this help message and exit

Base arguments:
  -ITS ITS    Sequence of ITS1 in fasta format.
  -type TYPE  Type of its region: its1, its2, concatenate.
  -out OUT    Get output of your file.

Algorithm Parameter:
  -K K        K nearest neighbors. Default k=10.
  -e E        Epsilon distance for nearest neighbor. 0 <= ε <= 1. Default ε=0.5.
```
## Reference
Will be add later ...
 
