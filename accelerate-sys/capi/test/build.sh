#!/bin/bash

 clang++ -std=c++17 -Wall -Wno-nullability-completeness main.cpp -I../src ../build/libaccelerate_sparse.a -framework Accelerate
