conda create --name careTCR python=3.11.11 numpy pandas matplotlib seaborn scikit-learn
conda activate careTCR
conda install -n careTCR ipykernel --update-deps --force-reinstall
conda install -c conda-forge tensorflow
conda install -c conda-forge notebook
