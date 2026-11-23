# Install XMAS 1.1.5 from source

The tool is displayed as **XMAS**, while its Python import name is
`chimerax.xmas` (lowercase). Version 1.1.4 explicitly used `chimerax.XMAS`;
this case-only change from older installations can cause import conflicts on
Windows. Version 1.1.5 standardizes the lowercase spelling.

1. Download the source ZIP for the branch/version being tested and extract it
   into a new writable folder, outside Program Files. Do not overlay an old
   source or build directory. The folder must contain pyproject.toml and src.
2. In the ChimeraX command line, change to that exact folder before building.
   ChimeraX 1.10's TOML builder copies extra files relative to its current
   working directory before invoking setuptools.

```text
cd "C:/path/to/fresh/XMAS_source"
devel install .
```

If rebuilding a previously used source folder, run `devel clean .` after the
`cd` and before `devel install .` to remove old build outputs.

3. Restart ChimeraX completely before opening XMAS.
4. Open **Tools → Structure Analysis → XMAS**. Check that Analyze has an
   **All possible K–K distances…** action and Export has an **AlphaLink2**
   output choice.

## Install the built wheel

The source build produces `dist/chimerax_xmas-1.1.5-py3-none-any.whl`.
To explicitly reinstall it, use its full path in ChimeraX:

```text
toolshed install "C:/path/to/XMAS_source/dist/chimerax_xmas-1.1.5-py3-none-any.whl" reinstall true
```

Restart after installing. Do not run system-Python pip to install into ChimeraX.

## Check module loading and metadata

In **Tools → General → Shell** (the Python shell, not the ChimeraX command line):

```python
import importlib
import chimerax
import importlib.metadata as md
print("ChimeraX package paths:", list(chimerax.__path__))
print("XMAS metadata:", [(d.version, str(d.locate_file(""))) for d in md.distributions() if d.metadata.get("Name", "").lower().replace("_", "-") == "chimerax-xmas"])
m = importlib.import_module("chimerax.xmas")
print("Loaded module:", m.__file__)
```

The imported module should come from the intended user installation. Multiple
metadata versions may obscure which bundle ChimeraX selects. A warning about
editable metadata is not conclusive proof of an editable installation: the
ChimeraX 1.10 scanner emits it whenever it finds no package paths in a bundle's
file list and no editable finder/path file. Source-folder egg-info is another
possible source. Change the working directory away from the source and restart
before collecting diagnostic output.

Do not delete broad site-packages directories or files belonging to other
bundles. Diagnose exact metadata and module paths before cleaning duplicates.
A remaining warning associated with another metadata entry does not by itself
prove the new module failed to load.

## Verification limits

The namespace change was checked against ChimeraX 1.10's bundle-builder naming
logic. Export and distance-calculation tests pass independently of ChimeraX.
A live Windows/ChimeraX launch is still required to confirm the installation
and GUI; a successful wheel build alone does not establish that the tool loads.
