#!/bin/bash

./build.sh

pushd ./test

./test.sh

popd
