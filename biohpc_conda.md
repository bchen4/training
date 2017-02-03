###To install software under home directory using conda

First of all:
```shell
module load python/2.7.x-anaconda
```

Conda has different channels (specified by -c), which are the software hubs. Important ones:
* Default
* r (Also default, but conda separated python and R channels)
* conda-forge
* bioconda (with many useful bioinfo tools)

Then create an conda enviornment by using:
```shell
conda create -n new_numpy <package list, delmited by space>
conda create -c r -n new_R r-essentials
```

You need to activate a conda evn to use it
```shell
source activate new_R
#deactivate current conda evn
source deactivate
```

When you are inside a conda evn, you can use conda to install the packages you want (suppose we are in new_R):
```shell
conda install bioconductor
```

To view already installed packages:
```shell
conda list
```



You can also use *pip* for to install under home directory by using --user :
```python
pip install --user twobitreader
```


## Create new env in specific folder
```shell
conda create -p /project/shared/bicf_workflow_ref/chipseq_bchen4/ -c r r-essentials
#Add channels
conda config --add channels conda-forge
conda config --add channels r
conda config --add channels bioconda
pip install --user twobitreader
conda install -c r r-xml
```

Install bioconductor in R console:
```R
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("DiffBind","ChIPseeker"))
```

