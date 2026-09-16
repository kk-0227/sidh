#!/bin/bash

# 厳格な警告（-Wall -Wextra）とC++17/20規格でコンパイル
g++ -std=c++17 -O2 -Wall -Wextra -Iinclude src/main.cpp -o build/main

# ビルド成功時のみ実行
if [ $? -eq 0 ]; then
    ./build/main
fi