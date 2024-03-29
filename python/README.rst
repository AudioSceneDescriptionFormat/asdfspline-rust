Python bindings for ASDF splines
================================

Requirements: Rust compiler + Cargo (https://rustup.rs/).

Install `cbindgen <https://crates.io/crates/cbindgen>`__::

    cargo install cbindgen

Compile and install Python module::

    python -m pip install -e .

Run tests (using `pytest <https://docs.pytest.org/>`__)::

    python -m pytest


Examples
^^^^^^^^

A few Jupyter notebooks are available in the ``examples/`` directory.
You can also play with them online at:
https://mybinder.org/v2/gh/AudioSceneDescriptionFormat/asdfspline-rust/master?labpath=python/examples
