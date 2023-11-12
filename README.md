# asdfspline

Rust implementation of the spline type that's used in the Audio Scene
Description Format (ASDF), see
<https://AudioSceneDescriptionFormat.readthedocs.io/>.


## Requirements

* Rust compiler, Cargo (<https://rustup.rs/>)

The required Rust packages (a.k.a. "crates") are listed in the file
`Cargo.toml`.


## Tests

```
cargo test --workspace
```

There are further tests (using Python) in the `python/` directory.


## API Documentation

Run `cargo doc --workspace` in the main directory to create the documentation.
The generated HTML documentation can be accessed via
`target/doc/asdfspline/index.html` and `target/doc/asdfspline_ffi/index.html`.


## Python Bindings

Have a look in the `python/` directory.


## License

MIT OR Apache-2.0
