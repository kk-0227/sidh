#include <iostream>
#include "fp.hpp"

int main(){
    Fp a(80);
    Fp b(30);
    (a/0).print();
    a.mul(a.inv()).print();
}