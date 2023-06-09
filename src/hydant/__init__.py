from __future__ import annotations

__docformat__ = "restructuredtext"

# dependency assumptions
_hard_dependencies = ('geopandas', 'networkx')
_missing_dependencies = []

for _dep in _hard_dependencies:
    try:
        __import__(_dep)
    except ImportError as _e:  # pragma: no cover
        _missing_dependencies.append(f"{_dep}: {_e}")

if _missing_dependencies:  # pragma: no cover
    raise ImportError(
        "Unable to import required dependencies:\n" + "\n".join(_missing_dependencies)
    )
del _hard_dependencies, _dep, _missing_dependencies


# docstring using numpy style
__doc__ = """
hydant: HYdrological ANalysis Tool for Python
=============================================

**hydant** is a simple Python package 

"""
