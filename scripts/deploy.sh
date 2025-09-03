#!/bin/bash
cargo test
cargo clippy
cargo fmt --check