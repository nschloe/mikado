#!/bin/sh

clang-tidy *pp -checks=* -- -I/usr/include/trilinos -I/usr/include/eigen3 -std=c++11
