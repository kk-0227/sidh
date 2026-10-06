#!/bin/bash

# 厳格な警告（-Wall -Wextra）とC++17/20規格でコンパイル
g++ -std=c++17 -Wall -Wextra -O2 -Iinclude src/isogeny.cpp src/curve.cpp src/fp.cpp src/main.cpp -o build/main


# ビルド成功時のみ実行
if [ $? -eq 0 ]; then
    ./build/main
fi