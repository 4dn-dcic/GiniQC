# README

GiniQC takes cool files as input. For the user's convenience, we provide two scripts to convert bedpe and ncc files to the cool format for compatibility with GiniQC. These can be found under the `utilities` directory in this repository. In order to use these scripts, the user will need to download [pairix](https://github.com/4dn-dcic/pairix).

bedpe2cool usage for a mouse dataset:
```
bash bedpe2cool.sh input.bedpe mm10
```
[More on the bedpe format.](https://bedtools.readthedocs.io/en/latest/content/general-usage.html)

ncc2cool usage for a mouse dataset:
```
bash ncc2cool.sh input.ncc mm10
```
[More on the ncc format.](https://github.com/tjs23/nuc_processing/wiki/NCC-data-format)

For conversion from pairs format, see [cooler](https://github.com/mirnylab/cooler).
For conversion from hic format, see [hic2cool](https://github.com/4dn-dcic/hic2cool).