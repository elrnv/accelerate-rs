#!/usr/bin/env bash

# Build the C API

mkdir -p build
pushd build

clang -c -Wall -Wno-nullability-completeness -o lib.o -I ../src ../src/api.c && ar rc ../libaccelerate_sparse.a lib.o

popd
