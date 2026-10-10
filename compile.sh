#!/bin/bash
set -e

CXXFLAGS="-std=c++17 -Wall -Wextra -O2 -Iinclude"
# main.cpp 以外の実装ファイル(ライブラリ部分)
LIB=$(ls src/*.cpp | grep -v '/main\.cpp$')

mkdir -p build
g++ $CXXFLAGS $LIB src/main.cpp -o build/main                                # デモ
g++ $CXXFLAGS $LIB tests/test_all.cpp -o build/test_all                      # テスト
g++ $CXXFLAGS $LIB tests/dump_shared_secrets.cpp -o build/dump_shared_secrets # Sage 照合用

./build/test_all