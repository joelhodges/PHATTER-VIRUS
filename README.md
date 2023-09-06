# PHATTER-VIRUS

## Setup for TACC users
To set up PHATTER-VIRUS on TACC please execute the following commands:

```
ssh username@ls6.tacc.utexas.edu
cd $HOME
ln -s $WORK work-ls6
module load python3
pip3 install astropy --user
pip3 install seaborn --user
pip3 install specutils --user
pip3 install scikit-learn --user
cd /work/NUMBER/NAME/ls6
git clone https://github.com/grzeimann/PHATTER-VIRUS.git
```
