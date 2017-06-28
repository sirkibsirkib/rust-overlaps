# Rust implementation for Kucherov's algorithm of an ASPOP solver

## Rust and Cargo

To get started, you need to ensure rust and cargo are installed on your system. They can be found at https://www.rust-lang.org/en-US/.
Next, you simply clone this repo, and open this directory in your terminal.

You can use Cargo to compile the source code for a project into a binary native to your system. Use the `--release` argument to ensure the code is compiled _with_ all compiler optimizations if you intend to use the system for realistic applications (It makes a massive difference!)
```
cargo build --release