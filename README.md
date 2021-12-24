# Dependencies

Download libigl
```
    git clone https://github.com/libigl/libigl.git
```

Place this repo next to the libigl folder. For example, if libigl is installed in ```~/foo/libigl```, then 
```
    git clone https://github.com/Coelancanth/ARAP_libigl.git ~/foo/ARAP_libigl
```
# Compilation

Compile libigl first (since we need imgui, we need issue '''cmake''' to download dependencies )
```
cd libigl/
mkdir build
cd build
cmake ../
make
```

Then compile this project using the standard CMake routine
``` 
cd ARAP_libigl/
mkdir build
cd build
cmake ../
make
```
