cmake_minimum_required(VERSION 3.24)
project(lab1 C)

set(CMAKE_C_STANDARD 23)

add_executable(lab1 lab1.c return_codes.h)
add_executable(main main.c)
add_executable(notmain notmain.c)