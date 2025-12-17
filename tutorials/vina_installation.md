# Installing Vina on Conda Environment

**Step 1**: Create a conda environment

```shell
conda create -n antimicrdpa python=3.12
```

**Step 2**: Install packages from ```requirements.txt```

```shell
pip install -r requirements.txt
```

**Step 3**: Install Autodock Vina on WSL

```shell
sudo apt install autodock-vina
```

**Step 4**: Install Autodock Vina Python library from GitHub Repo

```shell
git clone https://github.com/ccsb-scripps/AutoDock-Vina
```
```shell
cd AutoDock-Vina
```
```shell
git checkout meeko_wrap
```
```shell
cd build/python
```
```shell
python setup.py clean
```
```shell
python setup.py build
```
```shell
python3 setup.py install
```
```shell
pip install .
```

**Step 5**: Verify Vina is installed

```shell
ls ~/miniconda3/envs/antimicrdpa/lib/python3.12/site-packages/vina
```

**Step 6**: Remove AutodockVina Githhub folder

```shell
cd ../../..
``` 

Ensure you are on the antimicrDPA folder

```shell
rm -rf AutoDock-Vina/
```

END