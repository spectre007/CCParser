# CCParser
A general parsing framework for outputs of Quantum Chemistry codes.

## Dependencies
CCParser requires
* Python 3.x
* Python modules `re`, `inspect`, `os` and `importlib`

## Usage
Create a CCParser object specifying the output and software:
```python
import CCParser as ccp
QC = ccp.Parser("./calculation.out", software="Q-Chem")
```
Access the parsed quantities from the `.results` member variable. Since it's
generally possible to parse the same quantity several times (e.g. geometry optimization),
each result is stored in a list-like container. Access the individual quantity either
by index or use the `get_first()`/`get_last()` function, e.g.
```python
res = QC.results
print(res.scf_energy[2])
print(res.scf_energy.get_last())
```
