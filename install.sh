# get to your target folder and git the respository

git clone https://github.com/schneebergerlab/GameteBinning_prac.git

# install necessary tools

cd GameteBinning_prac

srcall=`echo src_*/`
for src in ${srcall}; do cd ./${src}; echo ./${src}; make; cd ..; done

# add the dev_bin to your environmental PATH variable.

cwp=`pwd`
export PATH=${cwp}/dev_bin:$PATH
