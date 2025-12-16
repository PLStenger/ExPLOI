# ExPLOI
ExPLOI

### Installing pipeline :

First, open your terminal. Then, run these two command lines :

    cd -place_in_your_local_computer
    git clone https://github.com/PLStenger/ExPLOI.git

### Update the pipeline in local by :

    git pull
    
### If necessary, install softwares by :   

    cd 99_softwares/
    conda install -c bioconda fastqc
    conda install -c bioconda trimmomatic
    conda install -c bioconda multiqc

For install QIIME2, please refer to http://qiime.org/install/install.html

### Know the number of CPU (threads) of your computer (here for MacOs) :   

    sysctl hw.ncpu
    > hw.ncpu: 4

### Run by :   

    time nohup bash pipeline.sh &> pipeline.out
