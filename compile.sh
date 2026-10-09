#!/bin/bash
set -e

# C++17、警告を厳しめにしてコンパイル
mkdir -p build
g++ -std=c++17 -Wall -Wextra -O2 -Iinclude src/*.cpp -o build/main

./build/main