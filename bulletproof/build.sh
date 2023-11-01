# Copyright 2023 licenser.author
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# 编译链接库

# ringct > crypto > common > epee > easylogging

cd common
gcc -c -fpic *.c *.cpp -I ./ -I ../epee/ -I ../easylogging++
gcc -shared *.o -o ../lib/libcommon.so
# commom > epee
# common > easylogging++

cd ..
cd crypto
gcc -c -fpic *.c *.cpp ./crypto_ops_builder/*.c -I ./crypto_ops_builder/ -I ./ -I ../epee/ -I ../
gcc -shared *.o -o ../lib/libcrypto.so
# crypto > epee
# crypto > common

cd ..
cd easylogging++
gcc -c -fpic *.cc -I .
gcc -shared *.o -o ../lib/libeasylogging.so
# easylogging = 0

cd ..
cd epee
gcc -c -fpic *.c *.cpp -I ./ -I ./storages/ -I ../easylogging++/
gcc -shared *.o -o ../lib/libepee.so
# epee > easylogging

cd ..
cd ringct
gcc -c -fpic *.c *.cc *.cpp -I ./ -I ../ -I ../epee -I ../easylogging++/
gcc -shared *.o -o ../lib/libringct.so
# ringct > epee
# ringct > easylogging
# ringct > common 
# ringct > crypto

# 编译程序
cd ..
g++ main.cpp -o ./build/main -I ./ -I ./epee/ -Wl,--rpath=./lib -L./lib -lringct -lcrypto -lcommon -lepee -leasylogging -L /usr/local/lib -lboost_thread -lboost_chrono

# 执行程序
cd ..
./build/main