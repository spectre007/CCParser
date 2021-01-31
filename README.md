# CCParser
A general parsing framework for outputs of Quantum Chemistry codes.

## Dependencies
CCParser requires
* Python 3.x
* Python modules `re`, `inspect`, `os` and `importlib`
* optional, but recommended: `pandas`

## Installation
Until the parser is converted to a proper Python package, install the module by
downloading the CCParser folder (not its contents) into a directory in your `$PYTHONPATH`.

If `$PYTHONPATH` is not set on your machine, create a folder intended for custom Python packages
and add it to the `$PYTHONPATH` in your `.bashrc` (or `.zshrc` or other).
```bash
cd ~; mkdir .my_python_packages
echo 'PYTHONPATH=$PYTHONPATH:$HOME/.my_python_packages' >> .bashrc
```

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
