# CCParser
A general parsing framework for outputs of Quantum Chemistry codes.

## Dependencies
CCParser requires
* Python 3.x
* Python modules `re`, `inspect`, `os` and `importlib`
* optional, but recommended: `pandas`

## Usage
Create a CCParser object specifying the output and software:
```python
import CCParser as ccp

# short form
qc1 = ccp.Parser("calculation.out")

# full list of options
qc2 = ccp.Parser("calculation.out", software="Q-Chem", to_console=True, to_file=False, log_file="CCParser.log", to_json=False, json_file="CCParser.json")
```
Access the parsed quantities from the `.results` member variable. Since it's
generally possible to parse the same quantity several times (e.g. geometry optimization),
each result is stored in a list-like container. Access the individual quantity either
by index or use the `get_first()`/`get_last()` function, e.g.
```python
res = qc1.results
print(res.scf_energy[2])
print(res.scf_energy.get_last())
```
