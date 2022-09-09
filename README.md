# pyLLE

Quite of a new version, v4.0.0 introduce a lot of new stuff including: 
- Complete recoding of the core with better codding and commenting which hopefully makes it easier to implement custom features 
- Modification of interfacing with parameters and results through method attributes instead of bulky dictionaries
- Stability of the half step Fourier method, allowing to use a soliton solution as an original state for the LLE 
- Allowing arbitrary number of driving pump, according to Taheri et al. and our paper on Nature Communication (Moille et al. )
- Julia compatibility with version 1.1 and above 

## How to Cite Us?

Please, if you use this package and it helps you with your research adn publication, cite us in your paper.
You can cite our paper published in the Journal of Research of National Institute of Standards and Technology available [here](https://doi.org/10.6028/jres.124.012). 
Not only it allows us to have a better idea of new things people are interested in and how to keep improving the solver, but it also help us building a community where every body could help maintaining the solver to suit better the needs of everybody. 

## Example 

Interactvive example can be found [here](https://mybinder.org/v2/gh/gregmoille/pyLLE/4.0.0?labpath=example%2FTemporalDualPump.ipynb)



> Moille G, Li Q, Lu X, Srinivasan K (2019) pyLLE: A Fast and User Friendly Lugiato-Lefever Equation Solver. J Res Natl Inst Stan 124:124012. https://doi.org/10.6028/jres.124.012

You can also use the bibtex entry: 
```latex
@article{moille_pyLLE,
      author = {Gregory Moille and Qing Li and Xiyuan Lu and Kartik Srinivasan},
      title = {pyLLE: a Fast and User Friendly Lugiato-Lefever Equation Solver},
      year = {2019},
      volume = {124},
      pages = {124012},
      month = {2019-05-24},
      journal = {Journal of Research of NIST},
       doi = {https://doi.org/10.6028/jres.124.012},
     }
```
