# kendrew

Visualizing molecules from Python.

![goodsell](media/goodsell.png)
![1e57](media/1e57.png)
![1e57](media/1e57.gif)

## Installation

Requirements:

* [VisPy](http://vispy.org/)
* [GlumPy](https://github.com/glumpy/glumpy)
* [BioPython](http://biopython.org/)
* [NumPy](http://www.numpy.org/)

Install them with pip:

```bash
pip install numpy vispy biopython PyOpenGL
pip install -U git+https://github.com/glumpy/glumpy@54a7eab7c08fb
```

## Running

```bash
python3 goodsell.py data/1e57.pdb --chain
```
