# asdfspline

Rust implementation of the spline type that's used in the Audio Scene
Description Format (ASDF), see
<https://AudioSceneDescriptionFormat.readthedocs.io/>.

## Requirements

* Rust compiler, Cargo (<https://rustup.rs/>)

The required Rust packages (a.k.a. "crates") are listed in the file
`Cargo.toml`.

## API Documentation

Run `cargo doc --all` in the main directory to create the documentation.
The generated HTML documentation can be accessed via
[target/doc/asdfspline/index.html](index.html) and
[target/doc/asdfspline_ffi/index.html](../asdfspline_ffi/index.html).

## Updating the C Header File

The file `ffi/asdfspline.h` was generated with
[cbindgen](https://crates.io/crates/cbindgen) (`cargo install cbindgen`).
After changes in the API functions, it can be updated with

```
cbindgen ffi -o ffi/asdfspline.h
```

## Updating `README.md`

Using [cargo-readme](https://github.com/livioribeiro/cargo-readme) (`cargo install cargo-readme`):

```
cargo readme -o README.md
```

## License

See `Cargo.toml`.
