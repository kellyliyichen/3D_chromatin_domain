# Regulatory roles of three-dimensional structures of chromatin domains
This page contains major analysis code and example data for the following manuscript:
[https://doi.org/10.1101/2022.07.22.501196](https://doi.org/10.1101/2022.07.22.501196)

<p align="center" width="100%">
    <img src=img/TAD.jpg width=70% height=70%>
</p>


## Requirement
- [Python](https://www.python.org) (>=3.8.5)
- [NumPy](https://numpy.org) (>=1.19.1)
- [pandas](https://pandas.pydata.org/) (>=1.1.1)
- [scikit-learn ](https://scikit-learn.org/stable/) (>=0.23.2)
- [SciPy](https://scipy.org/) (>=1.5.2)
- [Matplotlib](https://matplotlib.org/) (>=3.3.2)


## Installation
```
git clone https://github.com/kellyliyichen/3D_chromatin_domain.git
```


## Instruction
- To calculate intra-scDomain ratio using imaging data, follow the instruction [here] (scripts/sc_imaging/README.md)
- To calculate intra-TAD ratio using SPRITE data, follow the instruction [here](scripts/bulk_sprite/sprite/README.md)
- To calculate intra-TAD ratio using Hi-C data, follow the instruction [here]](scripts/bulk_sprite/hic/README.md)
